package org.broadinstitute.hellbender.tools.exome.allelefraction;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
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
final public class AlleleFractionSegmenter {
    private final List<AllelicCount> data;
    private final List<SimpleInterval> positions;
    private final AllelicPanelOfNormals allelicPoN;
    private final int N;
    private int K;
    private AlleleFractionHiddenMarkovModel model;
    private double concentration;
    private final double[] distances;   //distances[n] is the n to n+1 distance

    private static double NEGLIGIBLE_POSTERIOR_FOR_M_STEP = 0.01;

    //private static final Logger logger = LogManager.getLogger(AlleleFractionSegmenter.class);
    private static final double MINIMUM_MEMORY_LENGTH = 1;
    private static final double MAXIMUM_MEMORY_LENGTH = 1e10;
    private static final double DEFAULT_MEMORY_LENGTH = 5e6;

    // (unnormalized) vague gamma prior on concentration
    private static final Function<Double, Double> PRIOR_ON_CONCENTRATION = alpha -> alpha*Math.exp(-alpha);
    private static final double MINIMUM_CONCENTRATION = 1e-4;
    private static final double MAXIMUM_CONCENTRATION = 5;
    private static final int MAX_INTEGRATION_EVALUATIONS = 1000;
    private static final UnivariateIntegrator UNIVARIATE_INTEGRATOR = new SimpsonIntegrator(1e-3, 1e-3, 5, 20);

    private static final double CONVERGENCE_THRESHOLD = 0.01;
    private static final double MEMORY_LENGTH_CONVERGENCE_THRESHOLD = 1e4;

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
        data = Utils.nonNull(acc).getCounts();
        N = data.size();
        K = initialNumStates;
        this.allelicPoN = Utils.nonNull(allelicPoN);
        concentration = 1;

        // evenly spaced minor fractions from 1/2 to 1/(2K)
        final double[] evenlySpacedAlleleFractions = IntStream.range(0, K).mapToDouble(n ->  ((double) K - n) / (2*K)).toArray();
        final double[] flatWeights = Collections.nCopies(K, 1.0/K).stream().mapToDouble(x->x).toArray();
        model = new AlleleFractionHiddenMarkovModel(evenlySpacedAlleleFractions, flatWeights, DEFAULT_MEMORY_LENGTH,
                allelicPoN, AlleleFractionInitializer.INITIAL_PARAMETERS);
        positions = data.stream().map(AllelicCount::getInterval).collect(Collectors.toList());
        distances = IntStream.range(0, N - 1)
                .mapToDouble(n -> AlleleFractionHiddenMarkovModel.calculateDistance(positions.get(n), positions.get(n + 1)))
                .toArray();
        model = learnModel();
    }

    public List<ModeledSegment> findSegments() {
        List<Integer> states = ViterbiAlgorithm.apply(data, positions, model);
        List<ModeledSegment> result = new ArrayList<>();

        int beginningOfCurrentSegment = 0;
        int currentState = states.get(0);
        String currentContig = positions.get(0).getContig();
        for (int n = 0; n <= N; n++) {
            //if contig or state has switched, make previous segment and start new one
            if (n == N || currentState != states.get(n) || !currentContig.equals(positions.get(n).getContig())) {
                final int previousSegmentStart = positions.get(beginningOfCurrentSegment).getStart();
                final int previousSegmentEnd = positions.get(n-1).getEnd();
                final int previousSegmentSnpCount = n - beginningOfCurrentSegment;
                final double previousSegmentMinorFraction = model.getMinorAlleleFraction(currentState);
                final SimpleInterval interval = new SimpleInterval(currentContig, previousSegmentStart, previousSegmentEnd);
                result.add(new ModeledSegment(interval, previousSegmentSnpCount, previousSegmentMinorFraction));
                if (n < N) {
                    currentState = states.get(n);
                    currentContig = positions.get(n).getContig();
                    beginningOfCurrentSegment = n;
                }
            }
        }
        return result;
    }

    @VisibleForTesting
    protected AlleleFractionHiddenMarkovModel learnModel() {
        int iteration = 0;
        while (iteration < 15) {
            final double momentum = 1 + 1.5 / (1.0 + (double) iteration++/20.0);
            final double[] oldWeights = model.getWeights();
            final double[] oldMinorFractions = model.getMinorFractions();
            final double oldMemoryLength = model.getMemoryLength();
            performEMIteration(momentum);
            final double[] newWeights = model.getWeights();
            final double[] newMinorFractions = model.getMinorFractions();
            final double newMemoryLength = model.getMemoryLength();
            if (oldWeights.length != newWeights.length) {
                continue;
            }
            for (int n = 0; n < oldWeights.length; n++) {
                if (Math.abs(oldWeights[n] - newWeights[n]) > CONVERGENCE_THRESHOLD
                        || Math.abs(oldMinorFractions[n] - newMinorFractions[n]) > CONVERGENCE_THRESHOLD
                        || Math.abs(oldMemoryLength - newMemoryLength) > MEMORY_LENGTH_CONVERGENCE_THRESHOLD) {
                    continue;
                }
            }
            //break;  // if no continue was hit model has converged
        }
        return model;
    }

    // update the model and the concentration parameter with a single EM step
    private void performEMIteration(final double momentum) {
        final Dirichlet priorOnWeights = Dirichlet.symmetricDirichlet(K, concentration);
        final ExpectationStep expectationStep = new ExpectationStep();
        final double[] newMinorAlleleFractions = reestimateMinorAlleleFractions(expectationStep);
        final double[] emWeights = new Dirichlet(priorOnWeights, expectationStep.pseudocounts()).effectiveMultinomialWeights();
        //TODO: model.getWeights() is normalized; emWeights isn't. . .
        final double[] newWeights = accelerateWeights(model.getWeights(), emWeights, momentum);
        final double newMemoryLength = reestimateMemoryLength(expectationStep);
        final AllelicBiasParameters newParameters = reestimateParameters(expectationStep);
        model = new AlleleFractionHiddenMarkovModel(newMinorAlleleFractions, newWeights, newMemoryLength, allelicPoN, newParameters);
        pruneUnusedComponents();
        K = model.hiddenStates().size();
        concentration = variationalBayesConcentration(newWeights);
    }


    //TODO: detect candidate components to prune heuristically, but perform a likelihood test to decide whether to absorb thwem into
    // TODO: neighbor or delete.
    //filter out components that have low weight and are too close to another component -- these will
    // die out eventually in EM, but very slowly, so we hasten their demise for quicker convergence
    private void pruneUnusedComponents() {
        final double minorFractionDistanceThreshold = 0.01;
        final double lowWeight = 0.02;
        final double superLowWeight = 5e-4;

        //reverse order to preserve the convention that 0th component is f = 1/2
        final SortedMap<Double, Double> sortedFractionsAndWeights = new TreeMap<>(Collections.reverseOrder());
        IntStream.range(0, K).forEach(n -> sortedFractionsAndWeights.put(model.getMinorAlleleFraction(n), model.getWeight(n)));

        final double[] sortedFractions = sortedFractionsAndWeights.keySet().stream().mapToDouble(x->x).toArray();
        final double[] sortedWeights = sortedFractionsAndWeights.values().stream().mapToDouble(x->x).toArray();

        // look for low-weight components that also have less weight than their neighbors
        // never prune 0th component, which is neutral f = 1/2
        final Set<Integer> componentsToPrune = IntStream.range(1, K)
                .filter(n -> sortedWeights[n] < lowWeight)
                .filter(n -> sortedWeights[n] < sortedWeights[n-1])
                .filter(n -> n + 1 < K && sortedWeights[n] < sortedWeights[n+1])
                .filter(n -> Math.abs(sortedFractions[n] - sortedFractions[n - 1]) < minorFractionDistanceThreshold
                || n + 1 < K && Math.abs(sortedFractions[n] - sortedFractions[n + 1]) < minorFractionDistanceThreshold)
                .boxed().collect(Collectors.toSet());

        IntStream.range(1, K).filter(n -> sortedWeights[n] < superLowWeight).forEach(componentsToPrune::add);

        final double[] prunedWeights = IntStream.range(0, K)
                .filter(n -> !componentsToPrune.contains(n)).mapToDouble(n -> sortedWeights[n]).toArray();
        final double[] prunedFractions = IntStream.range(0, K)
                .filter(n -> !componentsToPrune.contains(n)).mapToDouble(n -> sortedFractions[n]).toArray();
        model = new AlleleFractionHiddenMarkovModel(prunedFractions, prunedWeights, model.getMemoryLength(), allelicPoN, model.getParameters());
    }

    private static double[] accelerateWeights(final double[] oldWeights, final double[] emWeights, final double acceleration) {
        return IntStream.range(0, emWeights.length)
                .mapToDouble(n -> Math.exp(Math.log(oldWeights[n]) + acceleration*Math.log(emWeights[n]/oldWeights[n]))).toArray();
    }

    /**
     * Compute the effective value of the Dirichlet concentration parameter, which defines the prior on weights in
     * subsequent iterations.  This value is the expectation of the concentration with respect to it mean-field
     * variational Bayes posterior distribution.  This mean-field comprises the prior on concentration, the
     * concentration-dependent Dirichlet distribution normalization constant, and the Dirichlet likelihood of the
     * effective weights.
     *
     * @param effectiveWeights the values of the weights to be plugged in to the mean-field distribution on the concentration
     *                         in a variational Bayes setting
     * @return  new estimate of the concentration
     */
    private static double variationalBayesConcentration(final double[] effectiveWeights) {
        final int K = effectiveWeights.length;
        final double geometricMeanOfEffectiveWeights = Math.exp(Arrays.stream(effectiveWeights).map(Math::log).average().getAsDouble());

        final Function<Double, Double> distribution = alpha -> PRIOR_ON_CONCENTRATION.apply(alpha)
                * Math.pow(geometricMeanOfEffectiveWeights, alpha)  //likelihood
                * Math.exp(Gamma.logGamma(alpha) - K * Gamma.logGamma(alpha/K));    //normalization constant

        return UNIVARIATE_INTEGRATOR.integrate(MAX_INTEGRATION_EVALUATIONS, alpha ->  alpha * distribution.apply(alpha), MINIMUM_CONCENTRATION, MAXIMUM_CONCENTRATION)
                / UNIVARIATE_INTEGRATOR.integrate(MAX_INTEGRATION_EVALUATIONS, distribution::apply, MINIMUM_CONCENTRATION, MAXIMUM_CONCENTRATION);
    }

    private double reestimateMemoryLength(final ExpectationStep eStep) {
        final double current = model.getMemoryLength();
        final Function<Double, Double> objective = D -> IntStream.range(0, distances.length)
                .mapToDouble(n -> eStep.pForget(n)*Math.log(1 - Math.exp(-distances[n]/D)) - (1 - eStep.pForget(n))*(distances[n]/D))
                .sum();
        return OptimizationUtils.quickArgmax(objective, MINIMUM_MEMORY_LENGTH, MAXIMUM_MEMORY_LENGTH, current);
    }

    private double[] reestimateMinorAlleleFractions(final ExpectationStep eStep) {
        return IntStream.range(0, K).mapToDouble(state -> {
            final Function<Double, Double> objective = f -> IntStream.range(0, data.size())
                    .mapToDouble(n -> {
                        final double eStepPosterior = eStep.pStateAtPosition(state, n);
                        return eStepPosterior < NEGLIGIBLE_POSTERIOR_FOR_M_STEP ? 0 : eStepPosterior * model.logEmissionProbability(data.get(n), f);
                    })
                    .sum();
            return state == 0 ? 0.5 : OptimizationUtils.quickArgmax(objective, 0, 0.5, model.getMinorAlleleFraction(state));
        }).toArray();
    }

    private AllelicBiasParameters reestimateParameters(final ExpectationStep eStep) {
        final Function<AllelicBiasParameters, Double> emissionLogLikelihood = params -> {
            double logLikelihood = 0.0;
            for (int position = 0; position < N; position++) {
                for (int state = 0; state < K; state++) {
                    final double eStepPosterior = eStep.pStateAtPosition(state, position);
                    logLikelihood += eStepPosterior < NEGLIGIBLE_POSTERIOR_FOR_M_STEP ? 0 :eStepPosterior
                            * AlleleFractionHiddenMarkovModel.logEmissionProbability(data.get(position), model.getMinorAlleleFraction(state), params, allelicPoN);
                }
            }
            return logLikelihood;
        };

        final Function<Double, Double> meanBiasObjective = mean -> emissionLogLikelihood.apply(model.getParameters().copyWithNewMeanBias(mean));
        final double newMeanBias = OptimizationUtils.quickArgmax(meanBiasObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_MEAN_BIAS, model.getParameters().getMeanBias());

        final Function<Double, Double> biasVarianceObjective = variance -> emissionLogLikelihood.apply(model.getParameters().copyWithNewBiasVariance(variance));
        final double newBiasVariance = OptimizationUtils.quickArgmax(biasVarianceObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_BIAS_VARIANCE, model.getParameters().getBiasVariance());

        final Function<Double, Double> outlierProbabilityObjective = pOutlier -> emissionLogLikelihood.apply(model.getParameters().copyWithNewOutlierProbability(pOutlier));
        final double newOutlierProbability = OptimizationUtils.quickArgmax(outlierProbabilityObjective, 0, AlleleFractionInitializer.MAX_REASONABLE_OUTLIER_PROBABILITY, model.getParameters().getOutlierProbability());

        return new AllelicBiasParameters(newMeanBias, newBiasVariance, newOutlierProbability);
    }

    /**
     * Stores the results of the expectation (E) step in which we run the forward-backward algorithm
     * to obtain posterior probabilities of 1) each hidden state at each het position; 2) each possible pair
     * of hidden states at consecutive het positions; 3) hidden state memory being lost at each het position;
     * and the expected number of 4) transitions to each state with memory loss, totalled over all het positions
     */
    final private class ExpectationStep {
        private final double[][] pStateByPosition = new double[K][N];

        //probability that memory was lost in the n to n+1 transition
        private final double[] pForget = new double[N-1];

        // pseudocounts for each hidden state
        final double[] stateCounts = new double[K];


        public ExpectationStep() {
            final ForwardBackwardAlgorithm.Result<AllelicCount, SimpleInterval, Integer> fbResult =
                    ForwardBackwardAlgorithm.apply(data, positions, model);

            IntStream.range(0, K).forEach(state ->
                    pStateByPosition[state] = IntStream.range(0, N).mapToDouble(n -> Math.exp(fbResult.logProbability(n, state))).toArray());

            for (int position = 0; position < N-1; position++) {
                for (int from = 0; from < K; from++) {
                    for (int to = 0; to < K; to++) {
                        // probability that fromState -> toState transition occurred going from position n to n + 1
                        final double pTransition = Math.exp(fbResult.logProbability(position, Arrays.asList(from, to)));
                        if (to != from ) {
                            stateCounts[to] += pTransition;
                            pForget[position] += pTransition;
                        } else {
                            final double weight = model.getWeight(to);
                            final double priorPForget = 1 - Math.exp(-distances[position] / model.getMemoryLength());
                            // Bayes' Rule gives the probability that the state was forgotten given that toState == fromState
                            final double pForgetAndTransition = pTransition * (priorPForget * weight / ((1-priorPForget) + priorPForget * weight));
                            stateCounts[to] += pForgetAndTransition;
                            pForget[position] += pForgetAndTransition;
                        }
                    }
                }
            }
        }

        public double pForget(final int position) { return pForget[position]; }
        public double pStateAtPosition(final int state, final int position) { return pStateByPosition[state][position]; }
        public double[] pseudocounts() { return stateCounts; }
    }
}
