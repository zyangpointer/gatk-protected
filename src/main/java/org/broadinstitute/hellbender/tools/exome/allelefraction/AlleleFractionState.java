package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.utils.mcmc.AbstractParameterizedState;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The state of the allele fraction model, containing: <p>
 *      1.  minor allele fractions for each segment <p>
 *      2.  a global outlier probability <p>
 *      3.  the mean allelic bias <p>
 *      4.  the rate (mean / variance) of the allelic bias <p>
 * <p>
 * See docs/CNVs/CNV-methods.pdf for details.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionState extends AbstractParameterizedState {
    @SuppressWarnings("serial")
    public static final class MinorFractions extends ArrayList<Double> {
        public MinorFractions() { super(); }
        public MinorFractions(final List<Double> other) {
            super(new ArrayList<>(other));
        }
    }

    public static final String MEAN_BIAS_NAME = "mu";
    public static final String BIAS_VARIANCE_NAME = "var";
    public static final String P_OUTLIER_NAME = "pi";
    public static final String MINOR_FRACTIONS_NAME = "f";

    @Override
    protected <S extends AbstractParameterizedState> S copy(final Class<S> stateClass) {
        return stateClass.cast(new AlleleFractionState(meanBias(), biasVariance(), outlierProbability(),
                new MinorFractions(minorFractions())));
    }

    public AlleleFractionState(final double meanBias, final double biasVariance, final double outlierProbability,
                               final MinorFractions minorFractions) {
        super(Arrays.asList(
                new Parameter<>(MEAN_BIAS_NAME, meanBias),
                new Parameter<>(BIAS_VARIANCE_NAME, biasVariance),
                new Parameter<>(P_OUTLIER_NAME, outlierProbability),
                new Parameter<>(MINOR_FRACTIONS_NAME, minorFractions)));
    }

    public double meanBias() {
        return get(MEAN_BIAS_NAME, Double.class);
    }

    public double biasVariance() {
        return get(BIAS_VARIANCE_NAME, Double.class);
    }

    public double outlierProbability() {
        return get(P_OUTLIER_NAME, Double.class);
    }

    public MinorFractions minorFractions() {
        return get(MINOR_FRACTIONS_NAME, MinorFractions.class);
    }

    public AllelicBiasParameters getParameters() {
        return new AllelicBiasParameters(meanBias(), biasVariance(), outlierProbability());
    }

    public double minorFractionInSegment(final int segment) {
        return get(MINOR_FRACTIONS_NAME, MinorFractions.class).get(segment);
    }
}
