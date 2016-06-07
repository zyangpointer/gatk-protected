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
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.hmm.ForwardBackwardAlgorithm;
import org.broadinstitute.hellbender.utils.hmm.ViterbiAlgorithm;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
final public class CopyRatioSegmenter extends ClusteringGenomicHMMSegmenter<Double> {
    private double logCoverageStandardDeviation;


    private static final double DEFAULT_INITIAL_STANDARD_DEVIATION = 0.1;
    private static final double NEUTRAL_LOG_2_COPY_RATIO = 0.0;
    private static final double MAX_REASONABLE_STANDARD_DEVIATION = 1.0;
    private static final double MIN_LOG_2_COPY_RATIO = -5.0;
    private static final double MAX_LOG_2_COPY_RATIO = 5.0;

    /**
     *
     * @param initialNumStates  A liberal estimate of the number of hidden minor allele fraction values to
     *                          include in the model.  Hidden states are pruned as the model is learned.
     */
    public CopyRatioSegmenter(final int initialNumStates, final List<SimpleInterval> positions, List<Double> data) {
        super(initialNumStates, positions, data);
    }

    /**
     * evenly-spaced log-2 copy ratios
     * @param K the initial number of hidden states
     */
    @Override
    protected void initializeHiddenStateValues(final int K) {
        hiddenStateValues = GATKProtectedMathUtils.createEvenlySpacedPoints(-2, 2, K);
        hiddenStateValues[NEUTRAL_VALUE_INDEX] = NEUTRAL_LOG_2_COPY_RATIO;
    }

    @Override
    protected void initializeAdditionalParameters() {
        logCoverageStandardDeviation = DEFAULT_INITIAL_STANDARD_DEVIATION;
    }

    @Override
    protected ClusteringGenomicHMM<Double> makeModel() {
        return new CopyRatioHiddenMarkovModel(hiddenStateValues, weights, memoryLength, logCoverageStandardDeviation);
    }

    @Override
    protected void relearnAdditionalParameters(final ExpectationStep eStep) {
        final Function<Double, Double> emissionLogLikelihood = sd -> {
            double logLikelihood = 0.0;
            for (int position = 0; position < positions.size(); position++) {
                for (int state = 0; state < weights.length; state++) {
                    final double eStepPosterior = eStep.pStateAtPosition(state, position);
                    logLikelihood += eStepPosterior < NEGLIGIBLE_POSTERIOR_FOR_M_STEP ? 0 :eStepPosterior
                            * CopyRatioHiddenMarkovModel.logEmissionProbability(data.get(position), hiddenStateValues[state], sd);
                }
            }
            return logLikelihood;
        };

        logCoverageStandardDeviation = OptimizationUtils.quickArgmax(emissionLogLikelihood, 0,
                MAX_REASONABLE_STANDARD_DEVIATION, logCoverageStandardDeviation);
    }

    @Override
    protected double minHiddenStateValue() { return MIN_LOG_2_COPY_RATIO; }

    @Override
    protected double maxHiddenStateValue() { return  MAX_LOG_2_COPY_RATIO; }
}
