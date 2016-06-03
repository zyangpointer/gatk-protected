package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.tools.exome.SegmentedModel;
import org.broadinstitute.hellbender.utils.mcmc.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Given a {@link org.broadinstitute.hellbender.tools.exome.SegmentedModel} and counts of alt and ref reads over a list of het sites,
 * infers the minor allele fraction of each segment.  For example, a segment
 * with (alt,ref) counts (10,90), (11,93), (88,12), (90,10) probably has a minor allele fraction
 * somewhere around 0.1.  The model takes into account allelic bias due to mapping etc. by learning
 * a global gamma distribution on allelic bias ratios.
 *<p>
 * We define the bias ratio of each het locus to be the expected ratio of
 * mapped ref reads to mapped alt reads given equal amounts of DNA (that is, given
 * a germline het).  The model learns a common gamma distribution:
 *      bias ratio ~ Gamma(alpha = mu^2/sigma^2, beta = mu/sigma^2)
 * where mu and sigma^2 are the global mean and variance of bias ratios, and
 * alpha, beta are the natural parameters of the gamma distribution.
 *</p>
 * <p>
 * Each segment has a minor allele fraction f, and for each het within the locus
 * the number of alt reads is drawn from a binomial distribution with total count
 * n = #alt reads + #ref reads and alt probability f/(f + (1-f)*bias ratio) if the
 * locus is alt minor and (1-f)/(1-f + f*bias ratio) if the locus is ref minor.
 *</p>
 * <p>
 * Conceptually, the model contains latent variables corresponding to the bias ratio
 * and indicators for alt minor/ref minor at each het locus.  However, we integrate them
 * out and the MCMC model below only contains the minor allele fractions and
 * the three hyperparameters of the model: the two parameters of the gamma distribution
 * along with the global outlier probability.
 *</p>
 * See docs/CNVs/CNV-methods.pdf for a thorough description of the model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionModeller {
    private final SegmentedModel segmentedModel;
    final AllelicPanelOfNormals allelicPON;
    private final ParameterizedModel<AlleleFractionState, AlleleFractionData> model;
    private final List<Double> meanBiasSamples = new ArrayList<>();
    private final List<Double> biasVarianceSamples = new ArrayList<>();
    private final List<Double> outlierProbabilitySamples = new ArrayList<>();
    private final List<AlleleFractionState.MinorFractions> minorFractionsSamples = new ArrayList<>();
    private final int numSegments;

    public AlleleFractionModeller(final SegmentedModel segmentedModel) {
        this(segmentedModel, AllelicPanelOfNormals.EMPTY_PON);
    }

    public AlleleFractionModeller(final SegmentedModel segmentedModel, final AllelicPanelOfNormals allelicPON) {
        this.segmentedModel = segmentedModel;
        this.allelicPON = allelicPON;
        final AlleleFractionData data = new AlleleFractionData(segmentedModel, allelicPON);
        numSegments = data.numSegments();
        final AlleleFractionState initialState = new AlleleFractionInitializer(data).getInitializedState();

        // Initialization got us to the mode of the likelihood
        // if we approximate conditionals as normal we can guess the width from the curvature at the mode

        final AllelicBiasParameters initialParameters = initialState.getParameters();
        final AlleleFractionState.MinorFractions initialMinorFractions = initialState.minorFractions();

        final double meanBiasInitialStepSize = estimateWidthAtMode(meanBias ->
                AlleleFractionLikelihoods.logLikelihood(initialParameters.copyWithNewMeanBias(meanBias), initialMinorFractions, data), initialParameters.getMeanBias());
        final double biasVarianceInitialStepSize = estimateWidthAtMode(biasVariance ->
                AlleleFractionLikelihoods.logLikelihood(initialParameters.copyWithNewBiasVariance(biasVariance), initialMinorFractions, data), initialParameters.getBiasVariance());
        final double outlierProbabilityInitialStepSize = estimateWidthAtMode(outlierProbability ->
                AlleleFractionLikelihoods.logLikelihood(initialParameters.copyWithNewOutlierProbability(outlierProbability), initialMinorFractions, data), initialParameters.getOutlierProbability());
        final List<Double> minorFractionsInitialStepSizes = IntStream.range(0, numSegments).mapToDouble(segment ->
                estimateWidthAtMode(f -> AlleleFractionLikelihoods.segmentLogLikelihood(initialParameters, f, data.countsInSegment(segment), allelicPON), initialMinorFractions.get(segment)))
                .boxed().collect(Collectors.toList());

        final Sampler<Double, AlleleFractionState, AlleleFractionData> meanBiasSampler =
                new AlleleFractionSamplers.MeanBiasSampler(initialState, meanBiasInitialStepSize);
        final Sampler<Double, AlleleFractionState, AlleleFractionData> biasVarianceSampler =
                new AlleleFractionSamplers.BiasVarianceSampler(initialState, biasVarianceInitialStepSize);
        final Sampler<Double, AlleleFractionState, AlleleFractionData> outlierProbabilitySampler =
                new AlleleFractionSamplers.OutlierProbabilitySampler(initialState, outlierProbabilityInitialStepSize);
        final Sampler<AlleleFractionState.MinorFractions, AlleleFractionState, AlleleFractionData> minorFractionsSampler =
                new AlleleFractionSamplers.MinorFractionsSampler(initialState, minorFractionsInitialStepSizes);

        model = new ParameterizedModel.GibbsBuilder<>(initialState, data, AlleleFractionState.class)
                .addParameterSampler(AlleleFractionState.MEAN_BIAS_NAME, meanBiasSampler, Double.class)
                .addParameterSampler(AlleleFractionState.BIAS_VARIANCE_NAME, biasVarianceSampler, Double.class)
                .addParameterSampler(AlleleFractionState.P_OUTLIER_NAME, outlierProbabilitySampler, Double.class)
                .addParameterSampler(AlleleFractionState.MINOR_FRACTIONS_NAME, minorFractionsSampler, AlleleFractionState.MinorFractions.class)
                .build();
    }

    /**
     * Adds {@code numSamples - numBurnIn} Markov-Chain Monte-Carlo samples of the parameter posteriors (generated using
     * Gibbs sampling) to the collections held internally.  The current {@link AlleleFractionState} held internally is used
     * to initialize the Markov Chain.
     * @param numSamples    total number of samples per posterior
     * @param numBurnIn     number of burn-in samples to discard
     */
    public void fitMCMC(final int numSamples, final int numBurnIn) {
        //run MCMC
        final GibbsSampler<AlleleFractionState, AlleleFractionData> gibbsSampler = new GibbsSampler<>(numSamples, model);
        gibbsSampler.runMCMC();

        //update posterior samples
        meanBiasSamples.addAll(gibbsSampler.getSamples(AlleleFractionState.MEAN_BIAS_NAME, Double.class, numBurnIn));
        biasVarianceSamples.addAll(gibbsSampler.getSamples(AlleleFractionState.BIAS_VARIANCE_NAME, Double.class, numBurnIn));
        outlierProbabilitySamples.addAll(gibbsSampler.getSamples(AlleleFractionState.P_OUTLIER_NAME, Double.class, numBurnIn));
        minorFractionsSamples.addAll(gibbsSampler.getSamples(AlleleFractionState.MINOR_FRACTIONS_NAME, AlleleFractionState.MinorFractions.class, numBurnIn));
    }

    public List<Double> getmeanBiasSamples() {
        return Collections.unmodifiableList(meanBiasSamples);
    }
    public List<Double> getBiasVarianceSamples() {
        return Collections.unmodifiableList(biasVarianceSamples);
    }
    public List<Double> getOutlierProbabilitySamples() { return Collections.unmodifiableList(outlierProbabilitySamples); }
    public List<AlleleFractionState.MinorFractions> getMinorFractionsSamples() { return Collections.unmodifiableList(minorFractionsSamples); }

    public List<List<Double>> getMinorFractionSamplesBySegment() {
        final List<List<Double>> result = new ArrayList<>();
        for (int segment = 0; segment < numSegments; segment++) {
            final List<Double> thisSegment = new ArrayList<>();
            for (final AlleleFractionState.MinorFractions sample : minorFractionsSamples) {
                thisSegment.add(sample.get(segment));
            }
            result.add(thisSegment);
        }
        return result;
    }

    /**
     * Returns a list of {@link PosteriorSummary} elements summarizing the minor-allele-fraction posterior for each segment.
     * Should only be called after {@link AlleleFractionModeller#fitMCMC(int, int)} has been called.
     * @param credibleIntervalAlpha credible-interval alpha, must be in (0, 1)
     * @param ctx                   {@link JavaSparkContext} used for mllib kernel density estimation
     * @return                      list of {@link PosteriorSummary} elements summarizing the
     *                              minor-allele-fraction posterior for each segment
     */
    public List<PosteriorSummary> getMinorAlleleFractionsPosteriorSummaries(final double credibleIntervalAlpha, final JavaSparkContext ctx) {
        final int numSegments = segmentedModel.getSegments().size();
        final List<PosteriorSummary> posteriorSummaries = new ArrayList<>(numSegments);
        for (int segment = 0; segment < numSegments; segment++) {
            final int j = segment;
            final List<Double> minorFractionSamples =
                    minorFractionsSamples.stream().map(s -> s.get(j)).collect(Collectors.toList());
            posteriorSummaries.add(PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(minorFractionSamples, credibleIntervalAlpha, ctx));
        }
        return posteriorSummaries;
    }

    //guess the width of a probability distribution given the position of its mode using a gaussian approximation
    private static double estimateWidthAtMode(final Function<Double, Double> logPDF, final double mode) {
        final double EPSILON = Math.min(1e-6, Math.abs(mode)/2);    //adjust scale is mode is very near zero
        final double DEFAULT = 1.0; //should never be needed; only used if mode is not a mode
        final double secondDerivative = (logPDF.apply(mode + EPSILON) - 2*logPDF.apply(mode) + logPDF.apply(mode - EPSILON))/(EPSILON*EPSILON);
        return secondDerivative < 0 ? Math.sqrt(-1.0 / secondDerivative) : DEFAULT;
    }
}
