package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.BiConsumer;

public class GenotypeGVCFsIntegrationTest extends CommandLineProgramTest {

    private static final List<String> NO_EXTRA_ARGS = Collections.emptyList();
    private final String basePairGVCF = "gvcf.basepairResolution.gvcf";

    public GenotypeGVCFsIntegrationTest() throws IOException {
    }

    @Test
    public void testSingleGvcf() throws IOException {
        final File expected = getTestFile("gvcf.basepairResolution.output.vcf");
        assertGenotypeGVCFsGenotypesMatchExpected(getTestFile(basePairGVCF), expected);
    }

    @Test
    public void testLarge() throws IOException {
        final File expected = getTestFile("CEUTrio.20.21.expected.vcf");
        assertGenotypeGVCFsGenotypesMatchExpected(new File("src/test/resources/large/CEUTrio.20.21.gatk3.4.gvcf"), expected);
        assertGenotypeGVCFSAnnotationsMatchExpected(new File("src/test/resources/large/CEUTrio.20.21.gatk3.4.gvcf"), expected, NO_EXTRA_ARGS);
    }

    private static <T>  void assertForEachElementInLists(final List<T> actual, final List<T> expected, final BiConsumer<T, T> assertion){
        Assert.assertEquals(actual.size(), expected.size());
        for(int i = 0; i < actual.size(); i++){
            assertion.accept(actual.get(i),expected.get(i));
        }
    }

    @DataProvider(name="gvcfsToGenotype")
    public Object[][] gvcfsToGenotype(){
        return new Object[][] {
                {basePairGVCF, "gvcf.basepairResolution.output.vcf", NO_EXTRA_ARGS}, //base pair level gvcf
                {"testUpdatePGT.gvcf", "testUpdatePGT.output.vcf", NO_EXTRA_ARGS},   //
                {"gvcfExample1.vcf", "gvcfExample1.vcf.expected.vcf", NO_EXTRA_ARGS},
                {"combined_genotype_gvcf_exception.vcf","combined_genotype_gvcf_exception.output.vcf", NO_EXTRA_ARGS},
                {"combined_genotype_gvcf_exception.nocall.vcf","combined_genotype_gvcf_exception.output.vcf", NO_EXTRA_ARGS},
                {basePairGVCF, "ndaTest.expected.vcf", Collections.singletonList("-nda")},
                {basePairGVCF, "maxAltAllelesTest.expected.vcf", Arrays.asList("--maxAltAlleles", "1")},
                {basePairGVCF, "standardConfTest.expected.vcf", Arrays.asList("-stand_call_conf", "300", "-stand_emit_conf", "100")},
                {"spanningDel.combined.g.vcf", "spanningDel.combined.g.vcf.expected.vcf", NO_EXTRA_ARGS},
                {"spanningDel.delOnly.g.vcf", "spanningDel.delOnly.g.vcf.expected.vcf", NO_EXTRA_ARGS}
        };
    }


    @Test(dataProvider = "gvcfsToGenotype")
    public void testGenotypes(String input, String expected, List<String> extraArgs) throws IOException {
        assertGenotypeGVCFsGenotypesMatchExpected(getTestFile(input), getTestFile(expected), extraArgs);
    }

    @Test(dataProvider = "gvcfsToGenotype")
    public void testAnnotations(String input, String expected, List<String> extraArgs) throws IOException {
        assertGenotypeGVCFSAnnotationsMatchExpected(getTestFile(input), getTestFile(expected), extraArgs);
    }

    @Test
    public void testUpdatePGT() throws IOException {
        assertGenotypeGVCFsGenotypesMatchExpected(getTestFile("testUpdatePGT.gvcf"), getTestFile("testUpdatePGT.output.vcf"));
    }

    private void assertGenotypeGVCFsGenotypesMatchExpected(File input, File expected) throws IOException {
        assertGenotypeGVCFsGenotypesMatchExpected(input, expected, NO_EXTRA_ARGS);
    }

    private void assertGenotypeGVCFSAnnotationsMatchExpected(File input, File expected, List<String> additionalArguments) throws IOException {
        final File output = createTempFile("genotypegvcf","vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                .addVCF(input)
                .addOutput(output);

        additionalArguments.forEach(args::add);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> expectedVC = getVariantContexts(expected);
        final List<VariantContext> actualVC = getVariantContexts(output);
        assertForEachElementInLists(actualVC, expectedVC, (a, e) -> BaseTest.assertVariantContextsAreEqual(a, e,
                Arrays.asList("FS", //TODO There's a bug in GATK 3 computing FS
                         "QD", //TODO QD has a cap value and anything that reaches that is randomized.  It's difficult to reproduce the same random numbers accross gatk3 -> 4
                        "InbreedingCoeff"))); //TODO Inbreeding calculation changed between 3.4 and now
    }

    private void assertGenotypeGVCFsGenotypesMatchExpected(File input, File expected, List<String> additionalArguments) throws IOException {
        final File output = createTempFile("genotypegvcf","vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addReference(new File(b37_reference_20_21))
                    .addVCF(input)
                    .addOutput(output);

        additionalArguments.forEach(args::add);

        Utils.resetRandomGenerator();
        runCommandLine(args);

        final List<VariantContext> expectedVC = getVariantContexts(expected);
        final List<VariantContext> actualVC = getVariantContexts(output);
        assertForEachElementInLists(actualVC, expectedVC, BaseTest::assertVariantContextsHaveSameGenotypes);
    }

//
//    @Test(enabled = true)
//    public void testUpdatePGTStrandAlleleCountsBySample() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                baseTestString(" -V " + privateTestDir + "testUpdatePGT.vcf -A StrandAlleleCountsBySample", b37KGReference),
//                1,
//                Arrays.asList("a96b79e7c3689c8d5506083cb6d27390"));
//        executeTest("testUpdatePGT, adding StrandAlleleCountsBySample annotation", spec);
//    }
//
//    @Test(enabled = true)
//    public void combineSingleSamplePipelineGVCF() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                baseTestString(" -V:sample1 " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
//                        " -V:sample2 " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
//                        " -V:sample3 " + privateTestDir + "combine.single.sample.pipeline.3.vcf" +
//                        " -L 20:10,000,000-20,000,000", b37KGReference),
//                1,
//                Arrays.asList("bf3c1982ab6ffee410cb6a1fff6e7105"));
//        executeTest("combineSingleSamplePipelineGVCF", spec);
//    }
//
//    @Test(enabled = true)
//    public void testTetraploidRun() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                baseTestString(" -V:sample1 " + privateTestDir + "tetraploid-gvcf-1.vcf" +
//                        " -V:sample2 " + privateTestDir + "tetraploid-gvcf-2.vcf" +
//                        " -V:sample3 " + privateTestDir + "tetraploid-gvcf-3.vcf" +
//                        " -L " + privateTestDir + "tetraploid-gvcfs.intervals", b37KGReference),
//                1,
//                Arrays.asList("47d454936dc1f17cf4c4f84f02841346"));
//        executeTest("combineSingleSamplePipelineGVCF", spec);
//    }
//
//    @Test(enabled= true)
//    public void testMixedPloidyRun() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                baseTestString(" -V:sample1 " + privateTestDir + "haploid-gvcf-1.vcf" +
//                        " -V:sample2 " + privateTestDir + "tetraploid-gvcf-2.vcf" +
//                        " -V:sample3 " + privateTestDir + "diploid-gvcf-3.vcf" +
//                        " -L " + privateTestDir + "tetraploid-gvcfs.intervals", b37KGReference),
//                1,
//                Arrays.asList("5d79ea9de8ada8520d01284cf0c9f720"));
//        executeTest("combineSingleSamplePipelineGVCF", spec);
//    }
//
//    @Test(enabled = true)
//    public void combineSingleSamplePipelineGVCF_includeNonVariants() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                baseTestString(" -V:sample1 " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
//                        " -V:sample2 " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
//                        " -V:sample3 " + privateTestDir + "combine.single.sample.pipeline.3.vcf" +
//                        " --includeNonVariantSites -L 20:10,030,000-10,033,000 -L 20:10,386,000-10,386,500", b37KGReference),
//                1,
//                Arrays.asList("d69b43cac448f45218e77308fc01e9e6"));
//        executeTest("combineSingleSamplePipelineGVCF_includeNonVariants", spec);
//    }
//
//    @Test(enabled = true)
//    public void combineSingleSamplePipelineGVCFHierarchical() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                baseTestString(" -V " + privateTestDir + "combine.single.sample.pipeline.combined.vcf" +
//                        " -V:sample1 " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
//                        " -V:sample2 " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
//                        " -V:sample3 " + privateTestDir + "combine.single.sample.pipeline.3.vcf" +
//                        " -L 20:10,000,000-20,000,000", b37KGReference),
//                1,
//                Arrays.asList("7c93d82758bfb6e7efec257ef8a46217"));
//        executeTest("combineSingleSamplePipelineGVCFHierarchical", spec);
//    }
//
//    @Test(enabled = true)
//    public void combineSingleSamplePipelineGVCF_addDbsnp() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                baseTestString(" -V:sample1 " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
//                        " -V:sample2 " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
//                        " -V:sample3 " + privateTestDir + "combine.single.sample.pipeline.3.vcf" +
//                        " -L 20:10,000,000-11,000,000 --dbsnp " + b37dbSNP132, b37KGReference),
//                1,
//                Arrays.asList("5b60a7a9575ea83407aa61123960a0cc"));
//        executeTest("combineSingleSamplePipelineGVCF_addDbsnp", spec);
//    }
//
    @Test
    public void testJustOneSample() throws IOException {
        assertGenotypeGVCFsGenotypesMatchExpected(getTestFile("gvcfExample1.vcf"), getTestFile("gvcfExample1.vcf.expected.vcf"));
    }
//
//    @Test(enabled = true)
//    public void testSamplesWithDifferentLs() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T GenotypeGVCFs --no_cmdline_in_header -L 1:69485-69791 -o %s -R " + b37KGReference +
//                        " -V " + privateTestDir + "gvcfExample1.vcf" +
//                        " -V " + privateTestDir + "gvcfExample2.vcf",
//                1,
//                Arrays.asList("8407cb9a1ab34e705e5a54a0d4146d84"));
//        executeTest("testSamplesWithDifferentLs", spec);
//    }
//
    @Test
    public void testNoPLsException() throws IOException {
        // Test with input files with (1) 0/0 and (2) ./.
        final File genotypeOutput = getTestFile("combined_genotype_gvcf_exception.output.vcf");
        assertGenotypeGVCFsGenotypesMatchExpected(getTestFile("combined_genotype_gvcf_exception.vcf"), genotypeOutput);
        assertGenotypeGVCFsGenotypesMatchExpected(getTestFile("combined_genotype_gvcf_exception.nocall.vcf"), genotypeOutput);

    }

    @Test
    public void testNDA() throws IOException {
        assertGenotypeGVCFsGenotypesMatchExpected(getTestFile(basePairGVCF), getTestFile("ndaTest.expected.vcf"), Collections.singletonList("-nda"));
    }

    @Test
    public void testMaxAltAlleles() throws IOException {
        assertGenotypeGVCFsGenotypesMatchExpected(getTestFile(basePairGVCF), getTestFile("maxAltAllelesTest.expected.vcf"), Arrays.asList("--maxAltAlleles", "1"));
    }

    @Test
    public void testStandardConf() throws IOException {
        assertGenotypeGVCFsGenotypesMatchExpected(getTestFile(basePairGVCF), getTestFile("standardConfTest.expected.vcf"), Arrays.asList("-stand_call_conf", "300", "-stand_emit_conf", "100"));
    }
//
//    @Test
//    public void testStrandAlleleCountsBySample() throws IOException {
//        //HaplotypeCaller creates gVCF
//        final String CEUTRIO_BAM = validationDataLocation + "CEUTrio.HiSeq.b37.chr20.10_11mb.bam";
//        final WalkerTestSpec specHaplotypeCaller = new WalkerTestSpec(
//                "-T HaplotypeCaller --disableDithering " +
//                        String.format("-R %s -I %s ", b37KGReference, CEUTRIO_BAM) +
//                        "--no_cmdline_in_header -o %s -L 20:10130000-10134800 " +
//                        "-ERC GVCF --sample_name NA12878 -variant_index_type LINEAR " +
//                        "-variant_index_parameter 128000 -A StrandAlleleCountsBySample",
//                1, Arrays.asList("")
//        );
//        specHaplotypeCaller.disableShadowBCF(); //TODO: Remove when BaseTest.assertAttributesEquals() works with SC
//        final File gVCF = executeTest("testStrandAlleleCountsBySampleHaplotypeCaller", specHaplotypeCaller).getFirst().get(0);
//        List<String> gVCFList = getAttributeValues(gVCF, new String("SAC"));
//
//        //Use gVCF from HaplotypeCaller
//        final WalkerTestSpec spec = new WalkerTestSpec(
//                baseTestString(" -V " + gVCF.getAbsolutePath(), b37KGReference),
//                1,
//                Arrays.asList(""));
//        final File outputVCF = executeTest("testStrandAlleleCountsBySample", spec).getFirst().get(0);
//        List<String> outputVCFList = getAttributeValues(outputVCF, new String("SAC"));
//
//        // All of the SAC values in the VCF were derived from the gVCF
//        Assert.assertTrue(gVCFList.containsAll(outputVCFList));
//    }
//
//    @Test
//    public void testUniquifiedSamples() throws IOException {
//        //two copies of 5 samples; will also test InbreedingCoeff calculation for uniquified samples
//        WalkerTestSpec spec = new WalkerTestSpec(
//                baseTestString(" -V:sample1 " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
//                        " -V:sample1B " + privateTestDir + "combine.single.sample.pipeline.1.vcf" +
//                        " -V:sample2 " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
//                        " -V:sample2B " + privateTestDir + "combine.single.sample.pipeline.2.vcf" +
//                        " -V:combined1 " + privateTestDir + "combine.single.sample.pipeline.combined.vcf" +
//                        " -V:combined2 " + privateTestDir + "combine.single.sample.pipeline.combined.vcf" +
//                        " --uniquifySamples", b37KGReference),
//                1,
//                Arrays.asList("9a472c4e101fff4892efb9255c5cd8b3"));
//        executeTest("testUniquifiedSamples", spec);
//
//    }
//
//    /**
//     * Returns a list of attribute values from a VCF file
//     *
//     * @param vcfFile VCF file
//     * @param attributeName attribute name
//     *
//     * @throws IOException if the file does not exist or can not be opened
//     *
//     * @return list of attribute values
//     */
//    private List<String> getAttributeValues(final File vcfFile, final String attributeName) throws IOException {
//        final VCFCodec codec = new VCFCodec();
//        final FileInputStream s = new FileInputStream(vcfFile);
//        final LineIterator lineIteratorVCF = codec.makeSourceFromStream(new PositionalBufferedStream(s));
//        codec.readHeader(lineIteratorVCF);
//
//        List<String> attributeValues = new ArrayList<String>();
//        while (lineIteratorVCF.hasNext()) {
//            final String line = lineIteratorVCF.next();
//            Assert.assertFalse(line == null);
//            final VariantContext vc = codec.decode(line);
//
//            for (final Genotype g : vc.getGenotypes()) {
//                if (g.hasExtendedAttribute(attributeNamsxtendedAttribute(attributeName));
//                }
//            }
//        }
//
//        return attributeValues;
//    }
//
//    /**
//     * Section to test spanning deletions
//     */
//    @Test
//    public void testSpanningDeletions() throws IOException {
//        final String gvcf1 = privateTestDir + "spanningDel.1.g.vcf";
//        final String gvcf2 = privateTestDir + "spanningDel.2.g.vcf";
//        final String gvcf3 = privateTestDir + "spanningDel.3.g.vcf";
//
//        // create the genotyped VCF to use as a basis for comparison against all of the combined versions
//        // case 0: GenotypeGVCFs(1.g.vcf, 2.g.vcf, 3.g.vcf)
//        final WalkerTestSpec genotypeBase = new WalkerTestSpec(
//                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + gvcf1 + " -V " + gvcf2 + " -V " + gvcf3,
//                1,
//                Arrays.asList(""));
//        genotypeBase.disableShadowBCF();
//        final File genotypeBaseVCF = executeTest("genotypeBase", genotypeBase).getFirst().get(0);
//        final List<VariantContext> BASE_VARIANT_CONTEXTS = getVariantContexts(genotypeBaseVCF);
//
//        // case 1: GenotypeGVCFs(CombineGVCFs(1.g.vcf, 2.g.vcf), 3.g.vcf)
//        final WalkerTestSpec combine12 = new WalkerTestSpec(
//                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + gvcf1 + " -V " + gvcf2,
//                1,
//                Arrays.asList(""));
//        combine12.disableShadowBCF();
//        final File combined_gVCF12 = executeTest("combine12", combine12).getFirst().get(0);
//        final WalkerTestSpec genotype12_3 = new WalkerTestSpec(
//                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + combined_gVCF12.getAbsolutePath() + " -V " + gvcf3,
//                1,
//                Arrays.asList(""));
//        genotype12_3.disableShadowBCF();
//        final File genotype12_3VCF = executeTest("genotype12_3", genotype12_3).getFirst().get(0);
//        final List<VariantContext> VARIANT_CONTEXTS_12_3 = getVariantContexts(genotype12_3VCF);
//        testVCsAreEqual(BASE_VARIANT_CONTEXTS, VARIANT_CONTEXTS_12_3);
//
//        // case 2: GenotypeGVCFs(CombineGVCFs(CombineGVCFs(1.g.vcf, 2.g.vcf), 3.g.vcf))
//        final WalkerTestSpec combine12then3 = new WalkerTestSpec(
//                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + combined_gVCF12 + " -V " + gvcf3,
//                1,
//                Arrays.asList(""));
//        combine12then3.disableShadowBCF();
//        final File combined_gVCF12then3 = executeTest("combined_gVCF12then3", combine12then3).getFirst().get(0);
//        final WalkerTestSpec genotype12then3 = new WalkerTestSpec(
//                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + combined_gVCF12then3.getAbsolutePath(),
//                1,
//                Arrays.asList(""));
//        genotype12then3.disableShadowBCF();
//        final File genotype12then3VCF = executeTest("genotype12then3", genotype12then3).getFirst().get(0);
//        final List<VariantContext> VARIANT_CONTEXTS_12then3 = getVariantContexts(genotype12then3VCF);
//        testVCsAreEqual(BASE_VARIANT_CONTEXTS, VARIANT_CONTEXTS_12then3);
//
//        // case 3: GenotypeGVCFs(CombineGVCFs(CombineGVCFs(1.g.vcf, 3.g.vcf), 2.g.vcf))
//        final WalkerTestSpec combine13 = new WalkerTestSpec(
//                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + gvcf1 + " -V " + gvcf3,
//                1,
//                Arrays.asList(""));
//        combine13.disableShadowBCF();
//        final File combined_gVCF13 = executeTest("combine13", combine13).getFirst().get(0);
//        final WalkerTestSpec combine13then2 = new WalkerTestSpec(
//                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + combined_gVCF13 + " -V " + gvcf2,
//                1,
//                Arrays.asList(""));
//        combine13then2.disableShadowBCF();
//        final File combined_gVCF13then2 = executeTest("combined_gVCF13then2", combine13then2).getFirst().get(0);
//        final WalkerTestSpec genotype13then2 = new WalkerTestSpec(
//                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + combined_gVCF13then2.getAbsolutePath(),
//                1,
//                Arrays.asList(""));
//        genotype13then2.disableShadowBCF();
//        final File genotype13then2VCF = executeTest("genotype13then2", genotype13then2).getFirst().get(0);
//        final List<VariantContext> VARIANT_CONTEXTS_13then2 = getVariantContexts(genotype13then2VCF);
//        testVCsAreEqual(BASE_VARIANT_CONTEXTS, VARIANT_CONTEXTS_13then2);
//
//        // case 4: GenotypeGVCFs(CombineGVCFs(1.g.vcf, 2.g.vcf, 3.g.vcf))
//        final WalkerTestSpec combine123 = new WalkerTestSpec(
//                "-T CombineGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + gvcf1 + " -V " + gvcf2 + " -V " + gvcf3,
//                1,
//                Arrays.asList(""));
//        combine123.disableShadowBCF();
//        final File combined_gVCF123 = executeTest("combine123", combine123).getFirst().get(0);
//        final WalkerTestSpec genotype123 = new WalkerTestSpec(
//                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + combined_gVCF123.getAbsolutePath(),
//                1,
//                Arrays.asList(""));
//        genotype123.disableShadowBCF();
//        final File genotype123VCF = executeTest("genotype123", genotype123).getFirst().get(0);
//        final List<VariantContext> VARIANT_CONTEXTS_123 = getVariantContexts(genotype123VCF);
//        testVCsAreEqual(BASE_VARIANT_CONTEXTS, VARIANT_CONTEXTS_123);
//    }
//
    /**
     * Returns a list of VariantContext records from a VCF file
     *
     * @param vcfFile VCF file
     *
     * @throws IOException if the file does not exist or can not be opened
     *
     * @return list of VariantContext records
     */
    private static List<VariantContext> getVariantContexts(final File vcfFile) throws IOException {
        final VCFCodec codec = new VCFCodec();
        final FileInputStream s = new FileInputStream(vcfFile);
        final LineIterator lineIteratorVCF = codec.makeSourceFromStream(new PositionalBufferedStream(s));
        codec.readHeader(lineIteratorVCF);

        final List<VariantContext> VCs = new ArrayList<>();
        while ( lineIteratorVCF.hasNext() ) {
            final String line = lineIteratorVCF.next();
            Assert.assertFalse(line == null);
            VCs.add(codec.decode(line));
        }

        return VCs;
    }

    private static void testVCsAreEqual(final List<VariantContext> VCs1, final List<VariantContext> VCs2) {
        Assert.assertEquals(VCs1.size(), VCs2.size(), "number of Variant Contexts");
        for ( int i = 0; i < VCs1.size(); i++ ) {
            final VariantContext vc1 = VCs1.get(i);
            final VariantContext vc2 = VCs2.get(i);
            Assert.assertEquals(vc1.toStringDecodeGenotypes(), vc2.toStringDecodeGenotypes());
        }
    }
//
//
//    private static final String simpleSpanningDeletionsMD5 = "e8616a396d40b4918ad30189856ceb01";
//
//    @Test(enabled = true)
//    public void testSpanningDeletionsMD5() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + privateTestDir + "spanningDel.1.g.vcf -V " + privateTestDir + "spanningDel.2.g.vcf",
//                1,
//                Arrays.asList(simpleSpanningDeletionsMD5));
//        spec.disableShadowBCF();
//        executeTest("testSpanningDeletionsMD5", spec);
//    }
//
    @Test
    public void testSpanningDeletionsFromCombinedGVCF() throws IOException {
        assertGenotypeGVCFsGenotypesMatchExpected(getTestFile("spanningDel.combined.g.vcf"), getTestFile("spanningDel.combined.g.vcf.expected.vcf"));
    }
//
//    @Test(enabled = true)
//    public void testMultipleSpanningDeletionsMD5() {
//        WalkerTestSpec spec = new WalkerTestSpec(
//                "-T GenotypeGVCFs --no_cmdline_in_header -o %s -R " + b37KGReference +
//                        " -V " + privateTestDir + "spanningDel.1.g.vcf -V " + privateTestDir + "spanningDel.2.g.vcf -V " + privateTestDir + "spanningDel.3.g.vcf",
//                1,
//                Arrays.asList("1c418229117bc8f148a69eda9c496309"));
//        spec.disableShadowBCF();
//        executeTest("testMultipleSpanningDeletionsMD5", spec);
//    }
//
    @Test
    public void testSpanningDeletionDoesNotGetGenotypedWithNoOtherAlleles() throws IOException {
        //output is an empty vcf
        assertGenotypeGVCFsGenotypesMatchExpected(getTestFile("spanningDel.delOnly.g.vcf"), getTestFile("spanningDel.delOnly.g.vcf.expected.vcf"));
    }
}
