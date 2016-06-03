package org.broadinstitute.hellbender.tools.exome.segmentation;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.exome.ModeledSegment;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AllelicBiasParameters;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.tools.exome.segmentation.AlleleFractionHiddenMarkovModel;
import org.broadinstitute.hellbender.tools.exome.segmentation.AlleleFractionSegmenter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * Created by davidben on 5/16/16.
 */
public class AlleleFractionSegmenterUnitTest {
    @Test
    public void testSegmentation() {
        final double[] trueWeights = new double[] {0.2, 0.5, 0.3};
        final double[] trueMinorAlleleFractions = new double[] {0.1, 0.3, 0.5};
        final double trueMemoryLength = 1e5;

        // the initial guess model starts with these parameters -- we're just testing learning of weights, minor fractions, and memory length
        final AllelicBiasParameters trueParams = new AllelicBiasParameters(1.0, 0.01, 0.01);

        AlleleFractionHiddenMarkovModel trueModel = new AlleleFractionHiddenMarkovModel(trueMinorAlleleFractions, trueWeights,
                trueMemoryLength, AllelicPanelOfNormals.EMPTY_PON, trueParams);

        //TODO: extract method and perhaps put in AFSimulatedData
        // randomly set positions
        final Random random = new Random(271);
        final int chainLength = 10000;
        final List<SimpleInterval> positions = new ArrayList<>();
        int position = 1;
        for (int n = 0; n < chainLength; n++) {
            position += random.nextInt((int) (trueMemoryLength/4));
            final SimpleInterval interval = new SimpleInterval("chr1", position, position);
            positions.add(interval);
        }

        final List<Integer> trueStates = trueModel.generateHiddenStateChain(positions);


        //make observed data

        //TODO: extract method
        //translate to ApacheCommons' parametrization of the gamma distribution
        final double gammaShape = trueParams.getMeanBias() * trueParams.getMeanBias() / trueParams.getBiasVariance();
        final double gammaScale = trueParams.getBiasVariance() / trueParams.getMeanBias();
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(563));
        final GammaDistribution biasGenerator = new GammaDistribution(rng, gammaShape, gammaScale);

        final AllelicCountCollection counts = new AllelicCountCollection();
        for (int n = 0; n < positions.size(); n++) {
            final int numReads = 100;
            final double minorFraction = trueMinorAlleleFractions[trueStates.get(n)];
            final double bias = biasGenerator.sample();

            //flip a coin to decide alt minor (alt fraction = minor fraction) or ref minor (alt fraction = 1 - minor fraction)
            final double altFraction =  random.nextDouble() < 0.5 ? minorFraction : 1 - minorFraction;

            //the probability of an alt read is the alt fraction modified by the bias or, in the case of an outlier, random
            final double pAlt = random.nextDouble() < trueParams.getOutlierProbability() ? random.nextDouble()
                    : altFraction / (altFraction + (1 - altFraction) * bias);

            final int numAltReads = new BinomialDistribution(rng, numReads, pAlt).sample();
            final int numRefReads = numReads - numAltReads;
            counts.add(new AllelicCount(positions.get(n), numAltReads, numRefReads));
        }

        final AlleleFractionSegmenter segmenter = new AlleleFractionSegmenter(10, counts, AllelicPanelOfNormals.EMPTY_PON);
        final AlleleFractionHiddenMarkovModel learnedModel = segmenter.learnModel();
        final List<ModeledSegment> segments = segmenter.findSegments();
        final double[] segmentMinorFractions = segments.stream()
                .flatMap(s -> Collections.nCopies((int) s.getTargetCount(), s.getSegmentMean()).stream())
                .mapToDouble(x->x).toArray();
        final double[] truthMinorFractions = trueStates.stream().mapToDouble(n -> trueModel.getMinorAlleleFraction(n)).toArray();

        final double averageMinorFractionError = IntStream.range(0, truthMinorFractions.length)
                .mapToDouble(n -> Math.abs(segmentMinorFractions[n] - truthMinorFractions[n]))
                .average().getAsDouble();

        Assert.assertEquals(averageMinorFractionError, 0, 0.01);
    }



}