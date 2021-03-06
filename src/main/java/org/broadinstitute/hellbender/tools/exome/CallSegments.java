package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Calls segments as amplified, deleted, or copy number neutral given files containing tangent-normalized
 * read counts by target and a list of segments
 *
 * @author David Benjamin
 */
@CommandLineProgramProperties(
        summary = "Call segments as amplified, deleted, or copy number neutral given files containing tangent-normalized" +
                " read counts by target and a list of segments",
        oneLineSummary = "Call segments as amplified, deleted, or copy number neutral",
        programGroup = CopyNumberProgramGroup.class
)
public final class CallSegments extends CommandLineProgram{

    @Argument(
            doc = "Tangent-normalized read counts input file.",
            shortName = ExomeStandardArgumentDefinitions.TARGET_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.TARGET_FILE_LONG_NAME,
            optional = false
    )
    protected File tangentNormalizedCoverageFile;

    @Argument(
            doc = "Segments files",
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            optional = false
    )
    protected File segmentsFile;

    @Argument(
            doc = "Called segments output",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false
    )
    protected File outFile;

    @Argument(
            doc = "(Advanced) Assume that the input seg file is using the legacy format (e.g. generated by python ReCapSeg)."
            + "  NOTE:  The output will be in the format used by this program -- i.e. no preservation of legacy field names, etc.",
            shortName = ExomeStandardArgumentDefinitions.LEGACY_SEG_FILE_SHORT_NAME,
            fullName = ExomeStandardArgumentDefinitions.LEGACY_SEG_FILE_LONG_NAME,
            optional = true
    )
    protected boolean isLegacyFormatSegFile = false;

    @Override
    protected Object doWork() {
        final ReadCountCollection tangentNormalizedCoverage;
        try {
            tangentNormalizedCoverage = ReadCountCollectionUtils.parse(tangentNormalizedCoverageFile);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(tangentNormalizedCoverageFile, e);
        }
        List<ModeledSegment> segments = isLegacyFormatSegFile ? SegmentUtils.readModeledSegmentsFromLegacySegmentFile(segmentsFile) :
               SegmentUtils.readModeledSegmentsFromSegmentFile(segmentsFile);

        ReCapSegCaller.makeCalls(tangentNormalizedCoverage, segments);

        final String sample = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(tangentNormalizedCoverageFile);
        SegmentUtils.writeModeledSegmentFile(outFile, segments, sample);
        return "SUCCESS";
    }
}
