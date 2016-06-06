package org.broadinstitute.hellbender.tools.exome.segmentation;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionInitializer;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AllelicBiasParameters;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.segmentation.AlleleFractionHiddenMarkovModel;
import org.broadinstitute.hellbender.utils.Dirichlet;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
final public class AlleleFractionSegmenter extends ClusteringGenomicHMMSegmenter<AllelicCount> {
    private AllelicPanelOfNormals allelicPoN;
    private AllelicBiasParameters biasParameters;

    /**
     * Initialize the segmenter with its data and panel of normals, giving equal weight to a set of evenly-spaced
     * hidden minor allele fraction values.
     *
     * @param initialNumStates  A liberal estimate of the number of hidden minor allele fraction values to
     *                          include in the model.  Hidden states are pruned as the model is learned.
     * @param acc               The {@link AllelicCountCollection} data attached to this segmenter
     * @param allelicPoN        The {@link AllelicPanelOfNormals} attached to this segmenter
     */
    public AlleleFractionSegmenter(final int initialNumStates, final AllelicCountCollection acc,
                                   final AllelicPanelOfNormals allelicPoN) {
        super(initialNumStates, acc.getCounts().stream().map(c -> c.getInterval()).collect(Collectors.toList()), acc.getCounts());
        this.allelicPoN = Utils.nonNull(allelicPoN);
    }

    /**
     * evenly-spaced minor allele fractions going from 1/2 to 0
     * @param K the initial number of hidden states
     */
    @Override
    protected void initializeHiddenStateValues(final int K) {
        hiddenStateValues = IntStream.range(0, K).mapToDouble(n ->  ((double) K - n) / (2*K)).toArray();
    }

    @Override
    protected void initializeAdditionalParameters() {
        biasParameters = new AllelicBiasParameters(1.0, 1e-3, 1e-2);
    }

    @Override
    protected ClusteringGenomicHMM<AllelicCount> makeModel() {
        return new AlleleFractionHiddenMarkovModel(hiddenStateValues, weights, memoryLength, allelicPoN, biasParameters);
    }

    @Override
    protected void relearnAdditionalParameters(final ExpectationStep eStep) {
        final Function<AllelicBiasParameters, Double> emissionLogLikelihood = params -> {
            double logLikelihood = 0.0;
            for (int position = 0; position < positions.size(); position++) {
                for (int state = 0; state < weights.length; state++) {
                    final double eStepPosterior = eStep.pStateAtPosition(state, position);
                    logLikelihood += eStepPosterior < NEGLIGIBLE_POSTERIOR_FOR_M_STEP ? 0 :eStepPosterior
                            * AlleleFractionHiddenMarkovModel.logEmissionProbability(data.get(position), hiddenStateValues[state], params, allelicPoN);
                }
            }
            return logLikelihood;
        };

        final Function<Double, Double> meanBiasObjective = mean -> emissionLogLikelihood.apply(biasParameters.copyWithNewMeanBias(mean));
        final double newMeanBias = OptimizationUtils.quickArgmax(meanBiasObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_MEAN_BIAS, biasParameters.getMeanBias());

        final Function<Double, Double> biasVarianceObjective = variance -> emissionLogLikelihood.apply(biasParameters.copyWithNewBiasVariance(variance));
        final double newBiasVariance = OptimizationUtils.quickArgmax(biasVarianceObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_BIAS_VARIANCE, biasParameters.getBiasVariance());

        final Function<Double, Double> outlierProbabilityObjective = pOutlier -> emissionLogLikelihood.apply(biasParameters.copyWithNewOutlierProbability(pOutlier));
        final double newOutlierProbability = OptimizationUtils.quickArgmax(outlierProbabilityObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_OUTLIER_PROBABILITY, biasParameters.getOutlierProbability());

        biasParameters = new AllelicBiasParameters(newMeanBias, newBiasVariance, newOutlierProbability);
    }

    @Override
    protected double minHiddenStateValue() { return 0.0; }

    @Override
    protected double maxHiddenStateValue() { return  0.5; }
}
