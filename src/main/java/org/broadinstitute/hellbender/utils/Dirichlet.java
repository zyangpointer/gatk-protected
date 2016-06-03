package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.stream.IntStream;

/**
 * The Dirichlet distribution is a distribution on multinomial distributions: if pi is a vector of positive multinomial weights
 * such that sum_i pi[i] = 1, the Dirichlet pdf is P(pi) = [prod_i Gamma(alpha[i]) / Gamma(sum_i alpha[i])] * prod_i pi[i]^(alpha[i] - 1)
 *
 * The vector alpha comprises the sufficient statistics for the Dirichlet distribution.
 *
 * Since the Dirichlet is the conjugate prior to the multinomial, if one has a Dirichlet prior with concentration alpha
 * and observes each category i n_i times (assuming categories are drawn from a multinomial distribution pi)
 * the posterior is alpha_i -> alpha_i + n_i
 *
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class Dirichlet {
    final double[] alpha;

    public Dirichlet(final double[] alpha) {
        Utils.nonNull(alpha);
        Arrays.stream(alpha).forEach(a -> ParamUtils.isPositive(a, "Dirichlet parameters must be positive."));
        if (alpha.length < 1) {
            throw new IllegalArgumentException("Dirichlet parameters must have at least one element");
        }
        this.alpha = Arrays.copyOf(alpha, alpha.length);
    }

    public Dirichlet(final Dirichlet prior, final double[] counts) {
        Utils.nonNull(counts);
        Arrays.stream(counts).forEach(c -> ParamUtils.isPositive(c, "Counts must be positive."));
        if (counts.length != prior.size()) {
            throw new IllegalArgumentException("Counts and prior must have same length.");
        }
        alpha = IntStream.range(0, prior.size()).mapToDouble(n -> prior.alpha[n] + counts[n]).toArray();
    }

    public static Dirichlet symmetricDirichlet(final int numStates, final double concentration) {
        return new Dirichlet(Collections.nCopies(numStates, concentration/numStates).stream().mapToDouble(x->x).toArray());
    }

    // in variational Bayes one often needs the effective point estimate of a multinomial distribution with a
    // Dirichlet prior.  This value is not the mode or mean of the Dirichlet but rather the exp of the expected log weights.
    // note that these effective weights do not add up to 1.  This is fine because in any probabilistic model scaling all weights
    // amounts to an arbitrary normalization constant, but it's important to keep in mind because some classes may expect
    // normalized weights.  In that case the calling code must normalize the weights.
    public double[] effectiveMultinomialWeights() {
        final double digammaOfSum = Gamma.digamma(MathUtils.sum(alpha));
        return Arrays.stream(alpha).map(a -> Math.exp(Gamma.digamma(a) - digammaOfSum)).toArray();
    }

    public int size() { return alpha.length; }
}
