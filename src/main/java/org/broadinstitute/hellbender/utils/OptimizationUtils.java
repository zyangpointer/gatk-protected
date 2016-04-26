package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;

import java.util.function.Function;

/**
 * Created by davidben on 4/27/16.
 */
public class OptimizationUtils {
    //TODO: put in stringent and lax optimizers and maybe put in an enum argument to argmax
    protected static final MaxEval DEFAULT_MAX_EVAL = new MaxEval(1000);
    private static final double DEFAULT_RELATIVE_TOLERANCE = 0.001;
    private static final double DEFAULT_ABSOLUTE_TOLERANCE = 0.001;
    protected static final BrentOptimizer DEFAULT_OPTIMIZER = new BrentOptimizer(DEFAULT_RELATIVE_TOLERANCE, DEFAULT_ABSOLUTE_TOLERANCE);

    protected static final MaxEval QUICK_MAX_EVAL = new MaxEval(100);
    private static final double QUICK_RELATIVE_TOLERANCE = 0.01;
    private static final double QUICK_ABSOLUTE_TOLERANCE = 0.01;
    protected static final BrentOptimizer QUICK_OPTIMIZER = new BrentOptimizer(QUICK_RELATIVE_TOLERANCE, QUICK_ABSOLUTE_TOLERANCE);

    private OptimizationUtils() { }

    public static double argmax(final Function<Double, Double> function, final double min, final double max, final double guess) {
        final UnivariateObjectiveFunction objective = new UnivariateObjectiveFunction(x -> function.apply(x));
        final SearchInterval interval= new SearchInterval(min, max, guess);
        return DEFAULT_OPTIMIZER.optimize(objective, GoalType.MAXIMIZE, interval, DEFAULT_MAX_EVAL).getPoint();
    }

    public static double quickArgmax(final Function<Double, Double> function, final double min, final double max, final double guess) {
        final UnivariateObjectiveFunction objective = new UnivariateObjectiveFunction(x -> function.apply(x));
        final SearchInterval interval= new SearchInterval(min, max, guess);
        return QUICK_OPTIMIZER.optimize(objective, GoalType.MAXIMIZE, interval, QUICK_MAX_EVAL).getPoint();
    }

}
