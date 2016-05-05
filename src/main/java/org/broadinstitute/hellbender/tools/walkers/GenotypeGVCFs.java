package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.MinimalGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.GeneralPloidyFailOverAFCalculatorProvider;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Perform joint genotyping on gVCF files produced by HaplotypeCaller
 *
 * <p>
 * GenotypeGVCFs merges gVCF records that were produced as part of the Best Practices workflow for variant discovery
 * (see Best Practices documentation for more details) using the '-ERC GVCF' or '-ERC BP_RESOLUTION' mode of the
 * HaplotypeCaller, or result from combining such gVCF files using CombineGVCFs. This tool will produce correct genotype
 * likelihoods, re-genotype the newly merged record, and then re-annotate it.</p>
 *
 * <h3>Input</h3>
 * <p>
 * One HaplotypeCaller gVCF to genotype
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined, genotyped VCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T GenotypeGVCFs \
 *   -R reference.fasta \
 *   -V sample1.g.vcf \
 *   -O output.vcf
 * </pre>
 *
 * <h3>Caveat</h3>
 * <p>Only gVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call gVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</p>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle any ploidy (or mix of ploidies) intelligently; there is no need to specify ploidy
 * for non-diploid organisms.</p>
 *
 */
@CommandLineProgramProperties(summary = "genotype a gvcf file to produce a vcf", oneLineSummary = "genotype a gvcf file", programGroup = VariantProgramGroup.class)
public class GenotypeGVCFs extends VariantWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written", optional=false)
    private File outputFile;

    @Argument(fullName="includeNonVariantSites", shortName="allSites", doc="Include loci found to be non-variant after genotyping", optional=true)
    private boolean includeNonVariants = false;

    @ArgumentCollection
    private GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    /**
     * Which annotations to recompute for the combined output VCF file.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to recompute.  The single value 'none' removes the default annotations", optional=true)
    private List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{"InbreedingCoeff", "FisherStrand", "QualByDepth", "ChromosomeCounts", "StrandOddsRatio"}));

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    private DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // the genotyping engine
    private GenotypingEngine<?> genotypingEngine;
    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    private VariantContextWriter vcfWriter;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        vcfWriter = GATKVariantContextUtils.createVCFWriter(outputFile, getReferenceDictionary(), false);
        
        final SampleList samples = new IndexedSampleList(getHeaderForVariants().getGenotypeSamples()); //todo should this be getSampleNamesInOrder?
        // create the genotyping engine
        genotypingEngine = new MinimalGenotypingEngine(createUAC(), samples, new GeneralPloidyFailOverAFCalculatorProvider(genotypeArgs));
        // create the annotation engine
        annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(Collections.emptyList(), annotationsToUse, Collections.<String>emptyList(), dbsnp.dbsnp, Collections.emptyList());

        // take care of the VCF headers
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(getHeaderForVariants().getMetaDataInInputOrder());
        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions());
        headerLines.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());

        // add headers for annotations added by this tool
        headerLines.add(new VCFSimpleHeaderLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG, GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_NAME_DEPRECATED, "Represents any possible spanning deletion allele at this location"));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        if ( dbsnp.dbsnp != null  )
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);

        final Set<String> sampleNameSet = samples.asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(headerLines, sampleNameSet);
        vcfWriter.writeHeader(vcfHeader);

        logger.info("Notice that the -ploidy parameter is ignored in " + getClass().getSimpleName() + " tool as this is automatically determined by the input variant files");
    }

    @Override
    public void apply(VariantContext vc, ReadsContext reads, ReferenceContext ref, FeatureContext features ) {
        ref.setWindow(10,10);
        final VariantContext mergedVC = ReferenceConfidenceVariantContextMerger.merge(Collections.singletonList(vc), vc, includeNonVariants ? ref.getBase() : null, true, false);
        final VariantContext regenotypedVC = regenotypeVC(mergedVC, ref, features);
        if (regenotypedVC != null) {
            vcfWriter.add(regenotypedVC);
        }
    }

    private static VariantContext removeNonRefSymbolicAllele(final VariantContext vc) {
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        final List<Allele> alleles = new ArrayList<>(vc.getAlleles());
        alleles.remove(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        builder.alleles(alleles);
        return builder.make();
    }


    /**
     * Re-genotype (and re-annotate) a combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
    private VariantContext regenotypeVC(final VariantContext originalVC, final ReferenceContext ref, final FeatureContext features) {
        Utils.nonNull(originalVC);

        VariantContext result = originalVC;

        // only re-genotype polymorphic sites
        if ( result.isVariant() ) {
            VariantContext regenotypedVC = genotypingEngine.calculateGenotypes(result, GenotypeLikelihoodsCalculationModel.SNP, null);
            if ( ! isProperlyPolymorphic(regenotypedVC) ) {
                if (!includeNonVariants)
                    return null;
            }
            else {
                regenotypedVC = GATKVariantContextUtils.reverseTrimAlleles(regenotypedVC);
                result = addGenotypingAnnotations(originalVC.getAttributes(), regenotypedVC);
            }
        }

        // if it turned monomorphic then we either need to ignore or fix such sites
        boolean createRefGTs = false;
        if ( result.isMonomorphicInSamples() ) {
            if ( !includeNonVariants)
                return null;
            createRefGTs = true;
        }

        // Re-annotate and fix/remove some of the original annotations.
        // Note that the order of these actions matters and is different for polymorphic and monomorphic sites.
        // For polymorphic sites we need to make sure e.g. the SB tag is sent to the annotation engine and then removed later.
        // For monomorphic sites we need to make sure e.g. the hom ref genotypes are created and only then are passed to the annotation engine.
        // We could theoretically make 2 passes to re-create the genotypes, but that gets extremely expensive with large sample sizes.
        if ( createRefGTs ) {
            result = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, true)).make();
            result = annotationEngine.annotateContext(result, features, ref, null, a -> true);
        } else {
            result = annotationEngine.annotateContext(result, features, ref, null, a -> true);
            result = new VariantContextBuilder(result).genotypes(cleanupGenotypeAnnotations(result, false)).make();
        }

        return result;
    }

    /**
     * Determines whether the provided VariantContext has real alternate alleles
     *
     * @param vc  the VariantContext to evaluate
     * @return true if it has proper alternate alleles, false otherwise
     */
    private boolean isProperlyPolymorphic(final VariantContext vc) {
        return ( vc != null && !vc.isSymbolic() );
    }

    /**
     * Add genotyping-based annotations to the new VC
     *
     * @param originalAttributes the non-null annotations from the original VC
     * @param newVC the new non-null VC
     * @return a non-null VC
     */
    private VariantContext addGenotypingAnnotations(final Map<String, Object> originalAttributes, final VariantContext newVC) {
        // we want to carry forward the attributes from the original VC but make sure to add the MLE-based annotations
        final Map<String, Object> attrs = new HashMap<>(originalAttributes);
        attrs.put(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        attrs.put(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, newVC.getAttribute(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        if (newVC.hasAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY))
            attrs.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY, newVC.getAttribute(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));

        return new VariantContextBuilder(newVC).attributes(attrs).make();
    }


    /**
     * Cleans up genotype-level annotations that need to be updated.
     * 1. move MIN_DP to DP if present
     * 2. propagate DP to AD if not present
     * 3. remove SB if present
     * 4. change the PGT value from "0|1" to "1|1" for homozygous variant genotypes
     * 5. move GQ to RGQ if the site is monomorphic
     *
     * @param VC            the VariantContext with the Genotypes to fix
     * @param createRefGTs  if true we will also create proper hom ref genotypes since we assume the site is monomorphic
     * @return a new set of Genotypes
     */
    private List<Genotype> cleanupGenotypeAnnotations(final VariantContext VC, final boolean createRefGTs) {
        final GenotypesContext oldGTs = VC.getGenotypes();
        final List<Genotype> recoveredGs = new ArrayList<>(oldGTs.size());
        for ( final Genotype oldGT : oldGTs ) {
            final Map<String, Object> attrs = new HashMap<>(oldGT.getExtendedAttributes());

            final GenotypeBuilder builder = new GenotypeBuilder(oldGT);
            int depth = oldGT.hasDP() ? oldGT.getDP() : 0;

            // move the MIN_DP to DP
            if ( oldGT.hasExtendedAttribute("MIN_DP") ) {
                depth = Integer.parseInt((String)oldGT.getAnyAttribute("MIN_DP"));
                builder.DP(depth);
                attrs.remove("MIN_DP");
            }

            // move the GQ to RGQ
            if ( createRefGTs && oldGT.hasGQ() ) {
                builder.noGQ();
                attrs.put(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, oldGT.getGQ());
            }

            // remove SB
            attrs.remove("SB");

            // update PGT for hom vars
            if ( oldGT.isHomVar() && oldGT.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) ) {
                attrs.put(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, "1|1");
            }

            // create AD if it's not there
            if ( !oldGT.hasAD() && VC.isVariant() ) {
                final int[] AD = new int[VC.getNAlleles()];
                AD[0] = depth;
                builder.AD(AD);
            }

            if ( createRefGTs ) {
                final int ploidy = oldGT.getPloidy();
                final List<Allele> refAlleles = Collections.nCopies(ploidy,VC.getReference());

                //keep 0 depth samples as no-call
                if (depth > 0) {
                    builder.alleles(refAlleles);
                }

                // also, the PLs are technically no longer usable
                builder.noPL();
            }

            recoveredGs.add(builder.noAttributes().attributes(attrs).make());
        }
        return recoveredGs;
    }

    /**
     * Creates a UnifiedArgumentCollection with appropriate values filled in from the arguments in this walker
     * @return a complete UnifiedArgumentCollection
     */
    private UnifiedArgumentCollection createUAC() {
        final UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = new GenotypeCalculationArgumentCollection(genotypeArgs);
        return uac;
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null) {
            vcfWriter.close();
        }
    }
}
