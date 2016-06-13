package org.broadinstitute.hellbender.tools.exome.allelefraction;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.special.Beta;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;

import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * The Allele Fraction Model (after marginalizing latent parameters as described in docs/CNVs/CNV-methods.pdf)
 * contains the following parameters:
 *      1.  minor allele fractions for each segment
 *      2.  a global outlier probability
 *      3.  the mean allelic bias
 *      4.  the rate (mean / variance) of the allelic bias
 *
 * Note that 3 and 4 are hyperparameters specifying a gamma distribution prior on allelic bias -- the latent variables
 * for bias at each het site have been marginalized but the hyperparameters have not.
 *
 * The Allele Fraction Model samples the distribution of these parameters using Markov chain Monte Carlo and in principle
 * an initialization step is not necessary.  However, in practice this initialization finds the mode of the posterior
 * distributions in only a few iterations, whereas sampling would require many more.  Thus we greatly reduce the
 * number of burn-in samples that we must discard.
 *
 * The initialization is straightforward: first we set the minor fractions to reasonable guesses based on alt and ref
 * counts, assuming no allelic bias.  Then we numerically maximize the likelihood with respect to each parameter until
 * the likelihood converges to a maximum.  In practice this is the unique global maximum.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionInitializer {
    @VisibleForTesting
    protected static final double INITIAL_OUTLIER_PROBABILITY = 0.01;
    protected static final double INITIAL_MEAN_BIAS = 1.0;
    protected static final double INITIAL_BIAS_VARIANCE = 0.1;   // this is an overestimate, but starting small makes it slow for
    // mean bias to escape a bad initial guess
    protected static final double LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD = 0.5;
    protected static final int MAX_ITERATIONS = 50;

    //define maxima of search intervals for maximum likelihood -- parameter values above these would be ridiculous
    private static final double MAX_REASONABLE_OUTLIER_PROBABILITY = 0.1;
    protected static final double MAX_REASONABLE_MEAN_BIAS = 5.0;
    protected static final double MAX_REASONABLE_BIAS_VARIANCE = 1.0;

    //the minor allele fraction of a segment must be less than one half by definition
    private static final double MAX_MINOR_ALLELE_FRACTION = 0.5;

    // constants for Brent optimization
    protected static final MaxEval BRENT_MAX_EVAL = new MaxEval(1000);
    private static final double RELATIVE_TOLERANCE = 0.001;
    private static final double ABSOLUTE_TOLERANCE = 0.001;
    protected static final BrentOptimizer OPTIMIZER = new BrentOptimizer(RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE);

    private AlleleFractionState state;
    private static final Logger logger = LogManager.getLogger(AlleleFractionInitializer.class);

    /**
     * Do the initialization
     * @param data data
     */
    public AlleleFractionInitializer(final AlleleFractionData data) {
        state = new AlleleFractionState(INITIAL_MEAN_BIAS, INITIAL_BIAS_VARIANCE, INITIAL_OUTLIER_PROBABILITY, initialMinorFractions(data));
        double previousIterationLogLikelihood;
        double nextIterationLogLikelihood = Double.NEGATIVE_INFINITY;
        logger.info(String.format("Initializing allele-fraction model.  Iterating until log likelihood converges to within %.3f.",
                LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD));
        int iteration = 1;
        do {
            previousIterationLogLikelihood = nextIterationLogLikelihood;
            state = new AlleleFractionState(estimateMeanBias(data), estimateBiasVariance(data),
                    estimateOutlierProbability(data), estimateMinorFractions(data));
            nextIterationLogLikelihood = AlleleFractionLikelihoods.logLikelihood(state, data);
            logger.info(String.format("Iteration %d, model log likelihood = %.3f.", iteration, nextIterationLogLikelihood));
            iteration++;
        } while (iteration < MAX_ITERATIONS &&
                nextIterationLogLikelihood - previousIterationLogLikelihood > LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD);
    }

    /**
     *
     * @return the initialized state of the Allele Fraction Model
     */
    public AlleleFractionState getInitializedState() { return state; }

    /**
     *  Initialize minor fractions assuming no allelic bias <p></p>
     *
     * We integrate over f to get posterior probabilities (responsibilities) of alt / ref minor
     * that is, responsibility of alt minor is int_{0 to 1/2} f^a (1-f)^r df
     *          responsibility of ref minor is int_{0 to 1/2} f^r (1-f)^a df
     * these are proportional to I(1/2, a + 1, r + 1) and I(1/2, r + 1, a + 1),
     * respectively, where I is the (incomplete) regularized Beta function.
     * By definition these likelihoods sum to 1, ie they are already normalized. <p></p>
     *
     * Finally, we set each minor fraction to the responsibility-weighted total count of
     * reads in minor allele divided by total reads, ignoring outliers.
     */
    private AlleleFractionState.MinorFractions initialMinorFractions(final AlleleFractionData data) {
        final int numSegments = data.getNumSegments();
        final AlleleFractionState.MinorFractions result = new AlleleFractionState.MinorFractions(numSegments);
        for (int segment = 0; segment < numSegments; segment++) {
            double responsibilityWeightedMinorAlleleReadCount = 0.0;
            double responsibilityWeightedTotalReadCount = 0.0;
            for (final AllelicCount count : data.getCountsInSegment(segment)) {
                final int a = count.getAltReadCount();
                final int r = count.getRefReadCount();
                double altMinorResponsibility;
                try {
                    altMinorResponsibility = Beta.regularizedBeta(0.5, a + 1, r + 1);
                } catch (final MaxCountExceededException e) {
                    altMinorResponsibility = a < r ? 1.0 : 0.0; //if the special function can't be computed, give an all-or-nothing responsibility
                }
                responsibilityWeightedMinorAlleleReadCount += altMinorResponsibility * a + (1 - altMinorResponsibility) * r;
                responsibilityWeightedTotalReadCount += a + r;
            }

            // we achieve a flat prior via a single pseudocount for minor and non-minor reads, hence the  +1 and +2
            result.add((responsibilityWeightedMinorAlleleReadCount + 1)/(responsibilityWeightedTotalReadCount + 2));
        }
        return result;
    }

    private double estimateOutlierProbability(final AlleleFractionData data) {
        final UnivariateObjectiveFunction objective = new UnivariateObjectiveFunction(outlierProbability ->
                AlleleFractionLikelihoods.logLikelihood(state.shallowCopyWithProposedOutlierProbability(outlierProbability), data));
        final SearchInterval searchInterval = new SearchInterval(0.0, MAX_REASONABLE_OUTLIER_PROBABILITY, state.outlierProbability());
        return OPTIMIZER.optimize(objective, GoalType.MAXIMIZE, searchInterval, BRENT_MAX_EVAL).getPoint();
    }

    private double estimateMeanBias(final AlleleFractionData data) {
        final UnivariateObjectiveFunction objective = new UnivariateObjectiveFunction(meanBias ->
                AlleleFractionLikelihoods.logLikelihood(state.shallowCopyWithProposedMeanBias(meanBias), data));
        final SearchInterval searchInterval = new SearchInterval(0.0, MAX_REASONABLE_MEAN_BIAS, state.meanBias());
        return OPTIMIZER.optimize(objective, GoalType.MAXIMIZE, searchInterval, BRENT_MAX_EVAL).getPoint();
    }

    private double estimateBiasVariance(final AlleleFractionData data) {
        final UnivariateObjectiveFunction objective = new UnivariateObjectiveFunction(biasVariance ->
                AlleleFractionLikelihoods.logLikelihood(state.shallowCopyWithProposedBiasVariance(biasVariance), data));
        final SearchInterval searchInterval = new SearchInterval(0.0, MAX_REASONABLE_BIAS_VARIANCE, state.biasVariance());
        return OPTIMIZER.optimize(objective, GoalType.MAXIMIZE, searchInterval, BRENT_MAX_EVAL).getPoint();
    }

    private double estimateMinorFraction(final int segment, final AlleleFractionData data) {
        final UnivariateObjectiveFunction objective = new UnivariateObjectiveFunction(minorFraction -> {
            final AlleleFractionState proposal =
                    new AlleleFractionState(state.meanBias(), state.biasVariance(), state.outlierProbability(), minorFraction); //single-segment state
            return AlleleFractionLikelihoods.segmentLogLikelihood(proposal, 0, data.getCountsInSegment(segment), data.getPON());
        });
        final SearchInterval searchInterval = new SearchInterval(0.0, MAX_MINOR_ALLELE_FRACTION, state.segmentMinorFraction(segment));
        return OPTIMIZER.optimize(objective, GoalType.MAXIMIZE, searchInterval, BRENT_MAX_EVAL).getPoint();
    }

    private AlleleFractionState.MinorFractions estimateMinorFractions(final AlleleFractionData data) {
        return new AlleleFractionState.MinorFractions(
                IntStream.range(0, data.getNumSegments())
                .mapToDouble(segment -> estimateMinorFraction(segment, data))
                .boxed().collect(Collectors.toList()));
    }
}
