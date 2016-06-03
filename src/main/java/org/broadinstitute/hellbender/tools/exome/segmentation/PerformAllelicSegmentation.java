package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.tools.exome.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.segmentation.AlleleFractionSegmenter;

import java.io.File;
import java.util.List;

/**
 * Created by davidben on 5/23/16.
 */
@CommandLineProgramProperties(
        summary = "Segment genomic data into regions of constant minor allele fraction.  Only supports one sample input.",
        oneLineSummary = "Segment genomic data into regions of constant minor allele fraction",
        programGroup = CopyNumberProgramGroup.class
)
public class PerformAllelicSegmentation extends CommandLineProgram {
    protected static final String NUM_STATES_LONG_NAME = "numberOfStates";
    protected static final String NUM_STATES_SHORT_NAME = "numStates";

    protected static final String SAMPLE_LONG_NAME = "sample";
    protected static final String SAMPLE_SHORT_NAME = "sample";


    @Argument(
            doc = "Input file for tumor-sample ref/alt read counts at normal-sample heterozygous-SNP sites (output of GetHetCoverage tool).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "Initial number of hidden allele fraction states",
            fullName = NUM_STATES_LONG_NAME,
            shortName = NUM_STATES_SHORT_NAME,
            optional = false
    )
    protected int initialNumStates;

    @Argument(
            doc = "Input file for allelic-bias panel of normals.",
            fullName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_SHORT_NAME,
            optional = true
    )
    protected File allelicPoNFile;

    @Argument(
            doc = "Output file for allele fraction segments.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File outputSegmentsFile;

    @Argument(
            doc = "Sample name.",
            fullName = SAMPLE_LONG_NAME,
            shortName = SAMPLE_SHORT_NAME,
            optional = false
    )
    protected String sampleName;

    @Override
    public Object doWork() {
        final AllelicPanelOfNormals allelicPoN =
                allelicPoNFile != null ? new AllelicPanelOfNormals(allelicPoNFile) : AllelicPanelOfNormals.EMPTY_PON;
        final AllelicCountCollection acc = new AllelicCountCollection(snpCountsFile);
        final AlleleFractionSegmenter segmenter = new AlleleFractionSegmenter(initialNumStates, acc, allelicPoN);
        segmenter.learnModel();
        final List<ModeledSegment> segments = segmenter.findSegments();
        SegmentUtils.writeModeledSegmentFile(outputSegmentsFile, segments, sampleName, true);

        return "SUCCESS";
    }
}
