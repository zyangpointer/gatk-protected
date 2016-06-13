package org.broadinstitute.hellbender.tools.exome.alleliccount;

import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountTableColumn.AllelicCountTableVerbosity;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;

import java.io.File;
import java.io.IOException;

/**
 * Writes {@link AllelicCountWithPhasePosteriors} instances to a tab-separated table file.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicCountWithPhasePosteriorsWriter extends AllelicCountWriter {
    public AllelicCountWithPhasePosteriorsWriter(final File file, final AllelicCountTableVerbosity verbosity) throws IOException {
        super(file, verbosity, PhasePosteriorsTableColumn.appendPhasePosteriorColumns(AllelicCountTableColumn.getColumns(verbosity)));
    }

    @Override
    protected void composeLine(final AllelicCount record, final DataLine dataLine) {
        super.composeLine(record, dataLine);
        Utils.validateArg(record.getClass().equals(AllelicCountWithPhasePosteriors.class), "Record must be of type AllelicCountWithPhasePosteriors.");
        composeLinePhasePosteriors((AllelicCountWithPhasePosteriors) record, dataLine);
    }

    /**
     * Compose the record for the phase posteriors.
     *
     * @param record the {@link AllelicCount} record
     * @param dataLine the {@link DataLine} to the composed
     */
    private static void composeLinePhasePosteriors(final AllelicCountWithPhasePosteriors record, final DataLine dataLine) {
        dataLine.append(formatProb(record.getRefMinorProb()))
                .append(formatProb(record.getAltMinorProb()))
                .append(formatProb(record.getOutlierProb()));
    }
}
