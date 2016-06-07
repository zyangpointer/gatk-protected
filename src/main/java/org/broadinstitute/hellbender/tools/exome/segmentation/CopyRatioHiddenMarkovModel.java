package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;
import java.util.stream.Collectors;

/**
 * The CopyRatio Hidden Markov Model is a generative model describing latent CNV states
 * with different copy ratios.
 *
 * Hidden states are essentially minor allele fractions f, but for convenience we represent them
 * as integers 0, 1, . . . K - 1 and store copy ratios in a corresponding array.
 *
 * The model contains a memory length parameter representing the prior probability for the minor
 * allele fraction state to be "forgotten" as a function of distance d between consecutive SNPs --
 * the probability to remember a state is exp(-d/memoryLength).  If the state is forgotten, a new state
 * is chosen with probabilities given by an array of weights.
 *
 * Thus our transition probabilities are P(i -> j) = exp(-d/D) delta_{ij} + (1 - exp(-d/D) weights[j]
 * where delta is the Kronecker delta and D is the memory length.
 *
 * Emission likelihoods as a function of copy ratio are assumed Gaussian with a shared global variance.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
final public class CopyRatioHiddenMarkovModel extends ClusteringGenomicHMM<Double> {
    private double logCoverageStandardDeviation;
    private List<AbstractRealDistribution> emissionDistributions;
    /**
     * @param log2CopyRatios array of log-2 copy ratios corresponding to the hidden states
     * @param weights array of (real-space, not log) prior probabilities of each hidden state
     *                when memory of the previous state is lost.  These may be unnormalized relative
     *                probabilities, which is useful when using variational Bayes.
     * @param memoryLength when consecutive SNPs are a distance d bases apart, the prior probability
     *                     for memory of the CNV state to be kept is exp(-d/memoryLength)
     */
    public CopyRatioHiddenMarkovModel(final double[] log2CopyRatios, final double[] weights,
                                      final double memoryLength, final double logCoverageStandardDeviation) {
        super(log2CopyRatios, weights, memoryLength);
        this.logCoverageStandardDeviation = logCoverageStandardDeviation;
        emissionDistributions = hiddenStates().stream()
                .map(n -> new NormalDistribution(log2CopyRatios[n], logCoverageStandardDeviation)).collect(Collectors.toList());
    }


    //TODO: make this robust to outliers with a uniform background or a Cauchy
    /**
     * Visible for the segmenter
     */
    public double logEmissionProbability(final Double data, final Integer state, final SimpleInterval position) {
        return emissionDistributions.get(state).logDensity(data);
    }

    public double logEmissionProbability(final Double data, final double copyRatio) {
        return logEmissionProbability(data, copyRatio, logCoverageStandardDeviation);
    }

    public static double logEmissionProbability(final Double data, final double copyRatio, final double standardDeviation) {
        return new NormalDistribution(null, copyRatio, standardDeviation).logDensity(data);
    }

    public double getCopyRatio(final int state) { return getHiddenStateValue(state); }
    public double[] getCopyRatios() { return getHiddenStateValues(); }
    public double getLogCoverageStandardDeviation() { return logCoverageStandardDeviation; }
}