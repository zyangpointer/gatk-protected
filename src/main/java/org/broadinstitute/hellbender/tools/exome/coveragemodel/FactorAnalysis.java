package org.broadinstitute.hellbender.tools.exome.coveragemodel;

import org.apache.commons.math3.linear.*;
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
    final int D;
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

    /**
     *
     * @param D dimension of observations
     * @param K dimension of latent space
     */
    public FactorAnalysis(final int D, final int K) {
        this.D = ParamUtils.isPositive(D, "Dimensionality of observations must be 1 or greater.");
        this.K = ParamUtils.isPositive(K, "Dimensionality of latent space must be 1 or greater.");
        ParamUtils.isPositiveOrZero(D-K, "Cannot have more latent dimensions than observed dimensions.");

        I = MatrixUtils.createRealIdentityMatrix(K);
        initializeParametersRandomly();

    }

    /**
     * learn parameters W, m, and Psi via expectation-maximization
     *
     * @param X NxD data matrix of observations
     * @param S NxD matrix of variance of observations -- variance of observation X_ij is S_ij
     */
    public void learn(final RealMatrix X, final RealMatrix S) {
        //TODO: check S is all positive

    }

    /**
     * One EM iteration including E step and M step
     * Note: this uses the approximations in which the variance of each sample is replaced by an average of sample variances
     * This is mainly for ease of implementation and we may use the exact M step later.
     */
    private void performOneExpectationMaximizationIteration() {

    }

    /*
     * Latent variables z are O(1) so each entry in the vector Wz is the inner product of K entries of W with K
     * entries of z.  If K is also O(1) this is a sum of K random-ish quantities that are all O(1), which has typical size
     * sqrt(K).  To get an O(1) result we thus need to initialize entries of W randomly with typical size 1/sqrt(K).
     *
     * We also need m and Psi to have O(1) entries for an O(1) result.
     */
    private void initializeParametersRandomly() {
        W = new Array2DRowRealMatrix(D, K);
        final double scale = 1/Math.sqrt(K);
        W.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) { return scale*rng.nextGaussian(); }
        });
        m = new ArrayRealVector(randomArray(D));
        Psi = new DiagonalMatrix(randomArray(D));
    }

    private double[] randomArray(final int size) {
        ParamUtils.isPositive(size, "Must have size > 0");
        return IntStream.range(0, D).mapToDouble(n -> rng.nextGaussian()).toArray();
    }

    private List<ExpectationStepResult> exactExpectationStepOfAllSamples(final RealMatrix X, final RealMatrix S) {
        return IntStream.range(0, X.getRowDimension())
                .mapToObj(s -> expectationStepForOneSample(X.getRowVector(s), S.getRowVector(s)))
                .collect(Collectors.toList());
    }

    private List<ExpectationStepResult> approximateExpectationStepOfAllSamples(final RealMatrix X, final RealMatrix S) {
        final int numSamples = S.getRowDimension();
        final RealVector averageSampleVariances = IntStream.range(0, numSamples)
                .mapToObj(s -> S.getRowVector(s))
                .reduce(new ArrayRealVector(S.getColumnDimension()), (a,b) -> a.add(b))
                .mapDivide(numSamples);
        final DiagonalMatrix totalVariance = Psi.add(new DiagonalMatrix(averageSampleVariances.toArray()));
    }

    public ExpectationStepResult expectationStepForOneSample(final RealVector observed, final RealVector variance) {
        if (observed.getDimension() != D) {
            throw new IllegalArgumentException(String.format("Input observed values must have length %d.", D));
        }

        final DiagonalMatrix sampleVariance = new DiagonalMatrix(variance.toArray());
        final DiagonalMatrix totalVariance = Psi.add(sampleVariance);
        final DiagonalMatrix precision = totalVariance.inverse();

        // G_i = (I + W^T (Psi+S_i)^(-1) W)^(-1)
        //this matrix is KxK so all operations are cheap.
        final RealMatrix G = new LUDecomposition(I.add(Wtranspose.multiply(precision).multiply(W))).getSolver().getInverse();
        return singleSampleEStepGivenPrecisionAndG(observed, precision, G);
    }


    /**
     * Calculate the E step for a single sample.  That is, calculate the 1st and 2nd moments of the Gaussian posterior
     * of the associated K-dimensional latent variable z
     * @param data D-dimensional data for a sample
     * @param precision precision matrix to be used in this calculation.  In exact mode this is
     *                  (Psi + Sigma_s)^-1, where Sigma_s is a diagonal matrix with jth entry equal to the
     *                  variance of the observed value data_j.  In approximate mode we use the average of Sigma_s over
     *                  all samples in order to have a shared precision and G matrix.
     * @param G         precomputed value of G = (I + W^T precision W)^(-1)
     * @return
     */
    private ExpectationStepResult singleSampleEStepGivenPrecisionAndG(final RealVector data,
                                                                      final DiagonalMatrix precision, final RealMatrix G) {
        // difference between data values and global mean m
        final RealVector residual = data.subtract(m);

        // E[z] = G W^T (Psi+S_i)^(-1) (x - m)
        final RealVector firstMoment = G.multiply(Wtranspose).multiply(precision).operate(residual);

        // E[z z^T] = G + E[z]E[z]^T
        final RealMatrix secondMoment = G.add(firstMoment.outerProduct(firstMoment));

        return new ExpectationStepResult(firstMoment, secondMoment);
    }

    /**
     * The E step of the EM optimization outputs the K-dimensional first moment vector E[z_i] and the KxK
     * second moment matrix E[z_i z_i^T] of each sample's latent vector with respect to its Gaussian posterior
     *
     * This class encapsulates these moments
     */
    private static final class ExpectationStepResult {
        public RealVector firstMoment;
        public RealMatrix secondMoment;

        public ExpectationStepResult(RealVector firstMoment, RealMatrix secondMoment) {
            this.firstMoment = firstMoment;
            this.secondMoment = secondMoment;
        }
    }
}
