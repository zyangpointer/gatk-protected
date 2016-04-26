package org.broadinstitute.hellbender.tools.exome.allelefraction;

/**
 * Encapsulates the global parameters of the allele fraction model: the mean and variance of the common prior on
 * allelic biases and the outlier probability.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class AllelicBiasParameters {
    private final double meanBias;
    private final double biasVariance;
    private final double outlierProbability;

    public AllelicBiasParameters(double meanBias, double biasVariance, double outlierProbability) {
        this.meanBias = meanBias;
        this.biasVariance = biasVariance;
        this.outlierProbability = outlierProbability;
    }

    public double getMeanBias() {
        return meanBias;
    }

    public double getBiasVariance() {
        return biasVariance;
    }

    public double getOutlierProbability() {
        return outlierProbability;
    }

    //get the gamma distribution alpha parameter
    public double getAlpha() {
        return meanBias * meanBias / biasVariance;
    }

    //get the gamma distribution beta parameter
    public double getBeta() {
        return meanBias / biasVariance;
    }

    public AllelicBiasParameters copyWithNewMeanBias(final double newMeanBias) {
        return new AllelicBiasParameters(newMeanBias, biasVariance, outlierProbability);
    }

    public AllelicBiasParameters copyWithNewBiasVariance(final double newBiasVariance) {
        return new AllelicBiasParameters(meanBias, newBiasVariance, outlierProbability);
    }

    public AllelicBiasParameters copyWithNewOutlierProbability(final double newOutlierProbability) {
        return new AllelicBiasParameters(meanBias, biasVariance, newOutlierProbability);
    }
}
