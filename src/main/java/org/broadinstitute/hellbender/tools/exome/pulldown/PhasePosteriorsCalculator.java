package org.broadinstitute.hellbender.tools.exome.pulldown;

import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionData;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionState;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountWithPhasePosteriorsCollection;

/**
 * Given an {@link AllelicCountCollection} and an {@link AlleleFractionState}
 * (representing the MAP model parameters for {@link AlleleFractionData} fit using {@link AlleleFractionModeller}),
 * calculates a new {@link AllelicCountWithPhasePosteriorsCollection} that also give
 * the probability for each het to be ref minor, alt minor, or an outlier, according to the model fit.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PhasePosteriorsCalculator {
    private final AllelicCountCollection counts;
    private final AlleleFractionState state;

    public PhasePosteriorsCalculator(final AllelicCountCollection counts, final AlleleFractionState state) {
        this.counts = counts;
        this.state = state;
    }

    public AllelicCountWithPhasePosteriorsCollection calculatePhasePosteriors() {
        return new AllelicCountWithPhasePosteriorsCollection();
    }
}
