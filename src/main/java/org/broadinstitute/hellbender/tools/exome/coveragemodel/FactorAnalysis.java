package org.broadinstitute.hellbender.tools.exome.coveragemodel;

import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Given N samples of D-dimensional data, represented as an NxD matrix X, where each measurement
 * X_ij is uncertain with variance S_ij, learn the following generative model:
 *
 * Each sample i is associated with a latent K-dimensional vector z_i (K << D) with isotropic Gaussian prior:
 * z_i ~ N(0, I)
 *
 * All samples share a common mean m and DxK matrix W that maps from K-dimensional latent space to the observed space.
 *
 * Measurements are taken from a normal distribution centered at the mapping from latent space:
 *
 * X_i ~ N(Wz_i + m, Psi)
 *
 * where Psi is a DxD diagonal matrix representing residual variance not explained by the latent variable z.
 *
 * We learn m, W, and Psi via the EM algorithm
 *
 * See CNV_methods.pdf for derivation of equations.
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class FactorAnalysis {
    int D;
    final int K;

    // the DxK matrix that maps from latent space to (variations in) observed space
    RealMatrix W;

    // cache the transpose
    RealMatrix Wtranspose;

    // the D-dimensional mean vector in observed space i.e. the expected value of observations when latent variables are zero
    RealVector m;

    // the residual variance, a DxD diagonal matrix
    DiagonalMatrix Psi;

    // the KxK identity matrix
    RealMatrix I;

    final static int RANDOM_SEED = 13;
    final Random rng = new Random(RANDOM_SEED);

    //TODO: javadoc
    public FactorAnalysis(final RealMatrix data, final RealMatrix variance, final int K) {
        if (data.getRowDimension() != variance.getRowDimension()) {
            throw new IllegalArgumentException(String.format("Data matrix has %d rows (samples) but variance matrix has %d.", data.getRowDimension(), variance.getRowDimension()));
        } else if (data.getColumnDimension() != variance.getColumnDimension()) {
            throw new IllegalArgumentException(String.format("Data matrix has %d columns (features) but variance matrix has %d.", data.getColumnDimension(), variance.getColumnDimension()));
        }
        this.K = ParamUtils.isPositive(K, "Must have positive number of features K.");
        I = MatrixUtils.createRealIdentityMatrix(K);
        initializeParameters(data);
        learn(data, variance);
    }

    /**
     * learn parameters W, m, and Psi via expectation-maximization
     *
     * @param data NxD data matrix of observations
     * @param variance NxD matrix of variance of observations -- variance of observation X_ij is S_ij
     */
    public void learn(final RealMatrix data, final RealMatrix variance) {
        //TODO: check S is all positive

    }

    private void performOneEMIteration(final RealMatrix data, final RealMatrix variance, final boolean approximate) {
        List<LatentVariablePosterior> latentVariablePosteriors = calculateLatentVariablePosteriors(data, variance, approximate);

    }

    /**
     * Initialize parameters with reasonable guesses.  Set m to the data mean.  Set Psi_jj to be the variance of the jth
     * data feature.  Initialize entries of W randomly with typical size 1/sqrt(K), which yields O(1) entries of Wz when
     * latent variables z are O(1).
     * @param data
     */
    private void initializeParameters(final RealMatrix data) {
        D = data.getColumnDimension();
        m = GATKProtectedMathUtils.averageOfRows(data);
        Psi = new DiagonalMatrix(GATKProtectedMathUtils.columnVariances(data));

        W = new Array2DRowRealMatrix(D, K);
        final double scale = 1/Math.sqrt(K);
        W.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) { return scale*rng.nextGaussian(); }
        });

    }

    private List<LatentVariablePosterior> calculateLatentVariablePosteriors(final RealMatrix data, final RealMatrix variance, final boolean approximate) {
        if (approximate) {
            final RealVector averageVariance = GATKProtectedMathUtils.averageOfRows(variance);
            Pair<DiagonalMatrix, RealMatrix> precisionAndG = calculatePrecisionAndG(averageVariance);
            return GATKProtectedMathUtils.rowStream(data)
                    .map(d -> calculateLatentVariablePosterior(d, precisionAndG))
                    .collect(Collectors.toList());
        } else {
            return IntStream.range(0, data.getRowDimension())
                    .mapToObj(s -> calculateLatentVariablePosterior(data.getRowVector(s), variance.getRowVector(s)))
                    .collect(Collectors.toList());
        }
    }

    private Pair<DiagonalMatrix, RealMatrix> calculatePrecisionAndG(final RealVector sampleVariance) {
        final DiagonalMatrix totalVariance = Psi.add(new DiagonalMatrix(sampleVariance.toArray()));
        final DiagonalMatrix precision = totalVariance.inverse();
        // G_i = (I + W^T (Psi+S_i)^(-1) W)^(-1)  -- this matrix is KxK so all operations are cheap.
        final RealMatrix G = new LUDecomposition(I.add(Wtranspose.multiply(precision).multiply(W))).getSolver().getInverse();
        return new Pair<>(precision, G);
    }

    public LatentVariablePosterior calculateLatentVariablePosterior(final RealVector data, final RealVector variance) {
        return calculateLatentVariablePosterior(data, calculatePrecisionAndG(variance));
    }

    /**
     * Calculate the E step for a single sample.  That is, calculate the 1st and 2nd moments of the Gaussian posterior
     * of the associated K-dimensional latent variable z
     * @param data D-dimensional data for a sample
     * @param precisionAndG precision matrix and G matrix to be used in this calculation.
     *                       In exact mode precision this is
     *                      (Psi + Sigma_s)^-1, where Sigma_s is a diagonal matrix with jth entry equal to the
     *                      variance of the observed value data_j.  In approximate mode we use the average of Sigma_s over
     *                      all samples in order to have a shared precision and G matrix.
     *                       G is precomputed value of G = (I + W^T precision W)^(-1)
     * @return
     */
    private LatentVariablePosterior calculateLatentVariablePosterior(final RealVector data, final Pair<DiagonalMatrix, RealMatrix> precisionAndG) {
        final DiagonalMatrix precision = precisionAndG.getFirst();
        final RealMatrix G = precisionAndG.getSecond();

        // difference between data values and global mean m
        final RealVector residual = data.subtract(m);

        // E[z] = G W^T (Psi+S_i)^(-1) (x - m) and E[z z^T] = G + E[z]E[z]^T
        final RealVector firstMoment = G.operate(Wtranspose.operate(precision.operate(residual)));
        final RealMatrix secondMoment = G.add(firstMoment.outerProduct(firstMoment));
        return new LatentVariablePosterior(firstMoment, secondMoment);
    }

    /**
     * The E step of the EM optimization outputs the K-dimensional first moment vector E[z_i] and the KxK
     * second moment matrix E[z_i z_i^T] of each sample's latent vector with respect to its Gaussian posterior
     *
     * This class encapsulates these moments
     */
    private static final class LatentVariablePosterior {
        public RealVector firstMoment;
        public RealMatrix secondMoment;

        public LatentVariablePosterior(RealVector firstMoment, RealMatrix secondMoment) {
            this.firstMoment = firstMoment;
            this.secondMoment = secondMoment;
        }
    }
}
