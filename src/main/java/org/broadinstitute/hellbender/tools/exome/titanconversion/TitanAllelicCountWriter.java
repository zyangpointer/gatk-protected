package org.broadinstitute.hellbender.tools.exome.titanconversion;

import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

public class TitanAllelicCountWriter extends TableWriter<AllelicCount> {

    public TitanAllelicCountWriter(final File file) throws IOException {
        super(file, new TableColumnCollection(TitanAllelicCountTableColumns.FULL_COLUMN_NAME_ARRAY));
    }

    @Override
    protected void composeLine(AllelicCount record, DataLine dataLine) {

        // Chr	Position	Ref	RefCount	Nref	NrefCount	NormQuality
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getEnd())
                .append(record.getRefNucleotide().name())
                .append(record.getRefReadCount())
                .append(record.getAltNucleotide().name())
                .append(record.getAltReadCount())
                .append("");
    }
}
