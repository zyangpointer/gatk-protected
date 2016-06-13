package org.broadinstitute.hellbender.tools.exome.pulldown;

import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionData;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionState;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountWithPhasePosteriors;

/**
 * Given an {@link AllelicCountCollection} and an {@link AlleleFractionState} (representing the MLE model parameters for
 * {@link AlleleFractionData} fit using {@link AlleleFractionModeller}), calculates a new {@link AllelicCountCollection}
 * of {@link AllelicCountWithPhasePosteriors} that also give the probability for each het to be ref minor, alt minor,
 * or an outlier, according to the model fit.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class PhasePosteriorsCalculator {
}
