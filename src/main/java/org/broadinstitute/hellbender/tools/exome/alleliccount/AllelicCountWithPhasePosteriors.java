package org.broadinstitute.hellbender.tools.exome.alleliccount;

import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn.AllelicCountTableVerbosity;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

/**
 * Represents an {@link AllelicCount} with posterior probabilities for each site being ref minor, alt minor,
 * or an outlier, according to a model fit by {@link AlleleFractionModeller}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicCountWithPhasePosteriors extends AllelicCount {

    /* these are extra metadata and can be null */
    private final double refMinorProb;
    private final double altMinorProb;
    private final double outlierProb;

    /**
     * Construct the {@link AllelicCountWithPhasePosteriors} object.  If unnormalized probabilities are passed,
     * they will be normalized.
     * @param count         {@link AllelicCount} of any verbosity (see {@link AllelicCountTableVerbosity})
     * @param refMinorProb  posterior probability of ref-minor phase (can be unnormalized)
     * @param altMinorProb  posterior probability of alt-minor phase (can be unnormalized)
     * @param outlierProb   posterior probability of outlier phase (can be unnormalized)
     */
    public AllelicCountWithPhasePosteriors(final AllelicCount count,
                                           final double refMinorProb, final double altMinorProb, final double outlierProb) {
        super(Utils.nonNull(count, "AllelicCount cannot be null."));
        ParamUtils.isPositiveOrZero(altMinorProb, "Can't construct AllelicCountWithPhasePosteriors with negative ref-minor probability.");
        ParamUtils.isPositiveOrZero(refMinorProb, "Can't construct AllelicCountWithPhasePosteriors with negative alt-minor probability.");
        ParamUtils.isPositiveOrZero(outlierProb, "Can't construct AllelicCountWithPhasePosteriors with negative outlier probability.");
        final double norm = refMinorProb + altMinorProb + outlierProb;
        this.refMinorProb = refMinorProb / norm;
        this.altMinorProb = altMinorProb / norm;
        this.outlierProb = outlierProb / norm;
    }

    public double getRefMinorProb() {
        return refMinorProb;
    }

    public double getAltMinorProb() {
        return altMinorProb;
    }

    public double getOutlierProb() {
        return outlierProb;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        if (!super.equals(o)) {
            return false;
        }

        final AllelicCountWithPhasePosteriors that = (AllelicCountWithPhasePosteriors) o;

        return Double.compare(that.refMinorProb, refMinorProb) == 0 &&
                Double.compare(that.altMinorProb, altMinorProb) == 0 &&
                Double.compare(that.outlierProb, outlierProb) == 0;

    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        long temp;
        temp = Double.doubleToLongBits(refMinorProb);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(altMinorProb);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(outlierProb);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
}
