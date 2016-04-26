package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;
import org.broadinstitute.hellbender.utils.segmenter.SegmenterUnitTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

import static org.testng.Assert.*;

/**
 * Created by davidben on 5/23/16.
 */
public class PerformAllelicSegmentationIntegrationTest extends CommandLineProgramTest {
    private static final String TOOLS_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";
    private static final File ALLELIC_COUNTS_FILE = new File(TOOLS_TEST_DIRECTORY, "snps-for-allelic-integration.tsv");

    @Test
    public void testCommandLine() throws IOException {
        final File snpFile = ALLELIC_COUNTS_FILE;
        final File outputSegmentFile = createTempFile("segments", ".seg");
        final String sampleName = "SAMPLE";
        final int initialNumStates = 10;
        final String[] arguments = {
                "-" + ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME, snpFile.getAbsolutePath(),
                "-" + PerformAllelicSegmentation.SAMPLE_SHORT_NAME, sampleName,
                "-" + PerformAllelicSegmentation.NUM_STATES_SHORT_NAME, Integer.toString(initialNumStates),
                "-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME, outputSegmentFile.getAbsolutePath(),
        };
        runCommandLine(arguments);
    }

}