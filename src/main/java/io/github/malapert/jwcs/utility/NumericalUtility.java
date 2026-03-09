/* 
 * Copyright (C) 2014-2022 Jean-Christophe Malapert
 *
 * This file is part of JWcs.
 * 
 * JWcs is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
package io.github.malapert.jwcs.utility;

import io.github.malapert.jwcs.proj.exception.MathematicalSolutionException;
import io.github.malapert.jwcs.proj.exception.JWcsError;
import java.text.DecimalFormat;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;
import org.apache.commons.math3.util.FastMath;

/**
 * NumericalUtility class.
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 */
public final class NumericalUtility {
    
    /**
     * Bisection algorithm.
     */
    private static final BisectionSolver solverBisection = new BisectionSolver(1e-15);
    
    /**
     * Laguerre algorithm.
     */
    private static final LaguerreSolver solverLaguerre = new LaguerreSolver(1e-15);

    /**
     * Double tolerance for numerical precision operations sets to 1e-12.
     */
    public static final double DOUBLE_TOLERANCE = 1e-12;

    /**
     * Half PI value.
     */
    public static final double HALF_PI = FastMath.PI * 0.5d;

    /**
     * Two Pi value.
     */
    public static final double TWO_PI = FastMath.PI * 2.0d;

    /**
     * Compares two doubles.
     *
     * @param val1 first double
     * @param val2 second double
     * @param precision precision for the comparison
     * @return True when <code>val1</code> and <code>val2</code> are equals.
     */
    public static boolean equalValues(final double val1, final double val2, final double precision) {
        return FastMath.abs(val2 - val1) <= precision;
    }

    /**
     * Projects ra/dec on cartesian reference system.
     *
     * @param pos position in the sky
     * @return the position in the cartesian reference system
     */
    public static double[] radec2xyz(final double[] pos) {
        final double[] xyz = new double[3];
        xyz[0] = FastMath.cos(pos[1]) * FastMath.cos(pos[0]);
        xyz[1] = FastMath.cos(pos[1]) * FastMath.sin(pos[0]);
        xyz[2] = FastMath.sin(pos[1]);
        return xyz;
    }

    /**
     * Norm a vector.
     *
     * @param pos position in the sky
     * @return the norm of the position.
     */
    public static double normVector(final double[] pos) {
        return FastMath.hypot(pos[0], pos[1]);
    }

    /**
     * Distance between two angles.
     *
     * @param pos1 first angle
     * @param pos2 second angle.
     * @return the distance between two angles
     */
    public static double distAngle(final double[] pos1, final double[] pos2) {
        final double[] xyzPos1 = radec2xyz(pos1);
        final double[] xyzPos2 = radec2xyz(pos2);
        double dot = xyzPos1[0] * xyzPos2[0] + xyzPos1[1] * xyzPos2[1] + xyzPos1[2] * xyzPos2[2];
        if (NumericalUtility.equalValues(dot, 0, 1e-13)) {
            dot = 0;
        }
        return NumericalUtility.aacos(dot / (normVector(xyzPos1) * normVector(xyzPos2)));
    }

    /**
     * Atan2 operation.
     *
     * @param n the ordinate coordinate
     * @param d the abscissa coordinate
     * @return the theta component of the point (r, theta) in polar coordinates
     * that corresponds to the point (x, y) in Cartesian coordinates.
     */
    public static double aatan2(final double n, final double d) {
        return aatan2(n, d, Double.NaN);
    }

    /**
     * Atan2 operation.
     *
     * @param n the ordinate coordinate
     * @param d the abscissa coordinate
     * @param defaultValue DefaultValue when atan2 is undefined
     * @return the theta component of the point (r, theta) in polar coordinates
     * that corresponds to the point (x, y) in Cartesian coordinates.
     */
    public static double aatan2(final double n, final double d, final double defaultValue) {
        return FastMath.abs(n) < DOUBLE_TOLERANCE && FastMath.abs(d) < DOUBLE_TOLERANCE ? defaultValue : FastMath.atan2(n, d);
    }

    /**
     * Asin operation.
     *
     * <p>Returns the arc sine of a value; the returned angle is in the range -pi/2
     * through pi/2. Special cases:
     * <ul>
     * <li>If the argument is NaN or its absolute value is greater than 1, then
     * the result is NaN.</li>
     * <li>If the argument is zero, then the result is a zero with the same sign
     * as the argument.</li>
     * </ul>
     *
     * @param v the value whose arc sine is returned
     * @return the arc sine of the argument.
     */
    public static double aasin(final double v) {
        if (equalValues(v, 1, DOUBLE_TOLERANCE)) {
            return FastMath.PI / 2;
        } else if (equalValues(v, -1, DOUBLE_TOLERANCE)) {
            return -FastMath.PI / 2;
        } else {
            return FastMath.asin(v);
        }
    }

    /**
     * Acos operation.
     *
     * @param v the value whose arc cosine is to be returned.
     * @return the value whose arc cosine is to be returned
     */
    public static double aacos(final double v) {
        if (FastMath.abs(v) > 1.) {
            return v < 0.0 ? FastMath.PI : 0.0;
        }
        return FastMath.acos(v);
    }

    /**
     * Normalizes the latitude.
     *
     * @param angle latitude in radians
     * @return the angle from -half_PI to half_PI
     */
    public static double normalizeLatitude(final double angle) {
        double resut = angle;
        if (Double.isInfinite(resut) || Double.isNaN(resut)) {
            throw new JWcsError("Infinite latitude");
        }
        if (equalValues(resut, HALF_PI)) {
            return HALF_PI;
        }
        if (equalValues(resut, -HALF_PI)) {
            return -HALF_PI;
        }

        if (resut > FastMath.PI) {
            resut -= TWO_PI;
        }

        if (resut < -FastMath.PI) {
            resut += TWO_PI;
        }
        return resut;
    }

    /**
     * normalizes longitude angle in radians.
     *
     * @param angle longitude in radians
     * @return the angle from 0 to 2PI
     */
    public static double normalizeLongitude(final double angle) {
        double result = angle;
        if (Double.isInfinite(angle) || Double.isNaN(angle)) {
            throw new JWcsError("Infinite longitude");
        }

        // avoid instable computations with very small numbers: if the
        // angle is very close to the graticule boundary, return +/-PI.
        // Bernhard Jenny, May 25 2010.
        if (FastMath.abs(result - 0) < DOUBLE_TOLERANCE) {
            return 0;
        }
        if (FastMath.abs(result - TWO_PI) < DOUBLE_TOLERANCE) {
            return TWO_PI;
        }

        while (result > TWO_PI) {
            result -= TWO_PI;
        }
        while (result < 0) {
            result += TWO_PI;
        }
        return result;
    }

    /**
     * Formats the number with a precision with 3 digits after the comma.
     *
     * @param number the number to format
     * @return the formatted number
     */
    public static String round(final double number) {
        final DecimalFormat df = new DecimalFormat("0.###");
        return df.format(number);
    }

    /**
     * Checks if the number is included in [min,max] with a numerical precision.
     *
     * @param number number to test
     * @param min minimum value
     * @param max maximum value
     * @param precision numerical precision
     * @return True when number is included in [min,max] otherwise False.
     */
    public static boolean isInInterval(final double number, final double min, final double max, final double precision) {
        return isInInterval(number, min, true, max, true, precision);
    }
    
    /**
     * Checks if the number is included in [min,max], [min,max[, ]min,max] or
     * ]min, max[ with a numerical precision.
     *
     * @param number number to test
     * @param min minimum value
     * @param minIsClosed the min value is included
     * @param max maximum value
     * @param maxIsClosed the max value is included
     * @param precision numerical precision
     * @return True when number is included in the interval otherwise False.
     */
    public static boolean isInInterval(final double number, final double min, final boolean minIsClosed, final double max, final boolean maxIsClosed, final double precision) {
        if (Double.isNaN(number)) {
            return false;
        }
        if (NumericalUtility.equalValues(number, min, precision)) {
            return minIsClosed;
        }
        if (NumericalUtility.equalValues(number, max, precision)) {
            return maxIsClosed;
        }
        return min < number && number < max;
    }    

    /**
     * Checks if the number is included in [min,max], [min,max[, ]min,max] or
     * ]min, max[ with a numerical precision of {@link NumericalUtility#DOUBLE_TOLERANCE}.
     *
     * @param number number to test
     * @param min minimum value
     * @param minIsClosed the min value is included
     * @param max maximum value
     * @param maxIsClosed the max value is included    
     * @return True when number is included in the interval otherwise False.
     */
    public static boolean isInInterval(final double number, final double min, final boolean minIsClosed, final double max, final boolean maxIsClosed) {
        if (Double.isNaN(number)) {
            return false;
        }
        if (NumericalUtility.equalValues(number, min)) {
            return minIsClosed;
        }
        if (NumericalUtility.equalValues(number, max)) {
            return maxIsClosed;
        }
        return min < number && number < max;
    }     
    
    /**
     * Compares two doubles with a numerical precision of
     * {@link NumericalUtility#DOUBLE_TOLERANCE}.
     *
     * @param val1 first double
     * @param val2 second double
     * @return True when <code>val1</code> and <code>val2</code> are equals.
     */
    public static boolean equalValues(final double val1, final double val2) {
        return equalValues(val1, val2, DOUBLE_TOLERANCE);
    }

    /**
     * Checks if the number is included in [min,max] with a numerical precision
     * of {@link NumericalUtility#DOUBLE_TOLERANCE}.
     *
     * @param number number to test
     * @param min minimum value
     * @param max maximum value
     * @return True when number is included in [min,max] otherwise False.
     */
    public static boolean isInInterval(final double number, final double min, final double max) {
        return isInInterval(number, min, max, DOUBLE_TOLERANCE);
    }

    /**
     * Calculates the matrix that represents a 3d rotation around the X axis.
     *
     * <p>Reference:<br>
     * ---------- <br>
     * Diebel, J. 2006, Stanford University, Representing Attitude: Euler
     * angles, Unit Quaternions and Rotation Vectors.
     * http://ai.stanford.edu/~diebel/attitude.html
     *
     * <p>Notes:<br>
     * ------<br>
     * Return the rotation matrix for a rotation around the X axis. This is a
     * rotation in the YZ plane. Note that we construct a new vector with: xnew
     * = R1.x In the literature, this rotation is usually called R1
     *
     * @param angle Rotation angle in degrees
     * @return A 3x3 matrix representing the rotation about angle around X axis.
     */
    public static RealMatrix rotX(final double angle) {
        final double angleRadians = FastMath.toRadians(angle);
        final double[][] array = {
            {1, 0, 0},
            {0, FastMath.cos(angleRadians), FastMath.sin(angleRadians)},
            {0, -FastMath.sin(angleRadians), FastMath.cos(angleRadians)}
        };
        return createRealMatrix(array);
    }

    /**
     * Calculates the matrix that represents a 3d rotation around the Y axis.
     *
     * <p>Reference:<br>
     * ----------<br>
     * Diebel, J. 2006, Stanford University, Representing Attitude: Euler
     * angles, Unit Quaternions and Rotation Vectors.
     * http://ai.stanford.edu/~diebel/attitude.html
     *
     * <p>Notes:<br>
     * ------<br>
     * Return the rotation matrix for a rotation around the X axis. This is a
     * rotation in the YZ plane. Note that we construct a new vector with: xnew
     * = R1.x In the literature, this rotation is usually called R1
     *
     * @param angle Rotation angle in degrees
     * @return A 3x3 matrix representing the rotation about angle around Y axis.
     */
    public static RealMatrix rotY(final double angle) {
        final double angleRadians = FastMath.toRadians(angle);
        final double[][] array = {
            {FastMath.cos(angleRadians), 0, -FastMath.sin(angleRadians)},
            {0, 1, 0},
            {FastMath.sin(angleRadians), 0, FastMath.cos(angleRadians)}
        };
        return createRealMatrix(array);
    }

    /**
     * Calculates the matrix that represents a 3d rotation around the Z axis.
     *
     * @param angle Rotation angle in degrees
     * @return A 3x3 matrix representing the rotation about angle around Z axis.
     */
    public static RealMatrix rotZ(final double angle) {
        final double angleRadians = FastMath.toRadians(angle);
        final double[][] array = {
            {FastMath.cos(angleRadians), FastMath.sin(angleRadians), 0},
            {-FastMath.sin(angleRadians), FastMath.cos(angleRadians), 0},
            {0, 0, 1}
        };
        return createRealMatrix(array);
    }

    /**
     * Creates an identity matrix.
     *
     * @param dimension matrix dimension
     * @return a matrix of dimension
     */
    public static RealMatrix createRealIdentityMatrix(final int dimension) {
        return MatrixUtils.createRealIdentityMatrix(3);
    }

    /**
     * A RealMatrix whose entries are the the values in the the input array.
     *
     * @param data input array
     * @return RealMatrix containing the values of the array
     */
    public static RealMatrix createRealMatrix(final double[][] data) {
        return MatrixUtils.createRealMatrix(data);
    }

    /**
     * Inverse matrix.
     *
     * @param matrix the matrix to inverse
     * @return the inverse matrix
     */
    public static RealMatrix inverse(final RealMatrix matrix) {
        return MatrixUtils.inverse(matrix);
    }

    /**
     * Computes the quadratic solution using the Laguerre's Method
     * for root finding of real coefficient polynomials.
     * 
     * <p>Once the roots are known, the nearest root from North pole is selected.
     * The first element of the coefficients array is the constant term. 
     * Higher degree coefficients follow in sequence. The degree of the 
     * resulting polynomial is the index of the last non-null element of the 
     * array, or 0 if all elements are null.
     * 
     * @param coefficients polynomial coefficients
     * @return the solution
     * @throws MathematicalSolutionException  No mathematical solution found
     * @see <a href="http://mathworld.wolfram.com/LaguerresMethod.html"> 
     * The Laguerre's method</a>
     */
    public static double computeQuatraticSolution(final double[] coefficients) throws MathematicalSolutionException {
        final Complex[] solutions = solverLaguerre.solveAllComplex(coefficients, 0);
        final Complex sol1 = solutions[0];
        final Complex sol2 = solutions[1];
        final double theta1 = NumericalUtility.aasin(sol1.getReal());
        final double theta2 = NumericalUtility.aasin(sol2.getReal());
        final boolean isTheta1Valid = NumericalUtility.isInInterval(theta1, -HALF_PI, HALF_PI);
        final boolean isTheta2Valid = NumericalUtility.isInInterval(theta2, -HALF_PI, HALF_PI);
        final double theta;
        if (isTheta1Valid && isTheta2Valid) {
            final double diffTheta1Pole = FastMath.abs(theta1 - HALF_PI);
            final double diffTheta2Pole = FastMath.abs(theta2 - HALF_PI);
            theta = diffTheta1Pole < diffTheta2Pole ? theta1 : theta2;
        } else if (isTheta1Valid) {
            theta = theta1;
        } else if (isTheta2Valid) {
            theta = theta2;
        } else {
            throw new MathematicalSolutionException("No mathematical solution found");
        }
        return theta;
    }

    /**
     * Solves for a zero root in the given interval for a given PolynomialFunction
     * object using a Bisection algorithm.
     * 
     * <p>A solver may require that the interval brackets a single zero root. 
     * Solvers that do require bracketing should be able to handle the case 
     * where one of the endpoints is itself a root.
     * 
     * @param maxEval maximum number of evaluations
     * @param f Polynomial function
     * @param min Lower bound for the interval
     * @param max Upper bound for the interval
     * @return a zero root
     * @see <a href="http://mathworld.wolfram.com/Bisection.html">Bisection algorithm</a>     
     */
    public static double computePolynomialSolution(final int maxEval, final Object f, final double min, final double max) {
        return solverBisection.solve(maxEval, (PolynomialFunction)f, min, max);
    }

    /**
     * Solves for a zero root in the given interval for a given function using
     * the bisection algorithm. 
     * 
     * <p>A solver may require that the interval brackets a single zero root. 
     * Solvers that do require bracketing should be able to handle the case 
     * where one of the endpoints is itself a root.
     * 
     * @param maxEval maximum number of evaluations
     * @param function function to solve
     * @param min Lower bound for the interval
     * @param max Upper bound for the interval
     * @return a zero root
     * @see <a href="http://mathworld.wolfram.com/Bisection.html">Bisection algorithm</a>
     */
    public static double computeFunctionSolution(final int maxEval, final UnivariateFunction function, final double min, final double max) {
        return solverBisection.solve(maxEval, function, min, max);
    }
  
    /**
     * Returns the degree of the polynomial function.
     * 
     * <p>f could be a polynomial function or an array representing the polynomial
     * coefficients.
     * 
     * @param f a polynomial function or an array representing the polynomial
     * coefficients
     * @return the degree of the polynomial function.
     * @throws JWcsError the type of object is not supported
     */
    public static int getPolynomialOrder(final Object f) {
        final int result;
        if (f instanceof PolynomialFunction fPolynomialFunction) {
            result = fPolynomialFunction.degree();
        } else if (f instanceof double[] fDouble) {
            final PolynomialFunction fNew = new PolynomialFunction(fDouble);
            result = fNew.degree();
        } else {
            throw new JWcsError("f is not a polynomialFunction or a an array of polynomial coefficients");
        }
        return result;
    }
    
    /**
     * Creates a polynomial function based on its coefficients.
     * 
     * <p>Constructs a polynomial with the given coefficients. The first element
     * of the coefficients array is the constant term. Higher degree coefficients
     * follow in sequence. The degree of the resulting polynomial is the index of
     * the last non-null element of the array, or 0 if all elements are null.
     * 
     * @param coefficients Polynomial coefficients
     * @return a polynomial function
     */
    public static Object createPolynomialFunction(final double[] coefficients) {
        return new PolynomialFunction(coefficients);
    }
    
    /**
     * Returns a copy of the coefficients array.
     * 
     * <p>Changes made to the returned copy will not affect the coefficients of 
     * the polynomial.
     * 
     * @param f the polynomial function of PolynomialFunction type
     * @return a fresh copy of the coefficients array
     */
    public static double[] getPolynomialCoefficients(final Object f) {
        final PolynomialFunction poly = (PolynomialFunction)f;
        return poly.getCoefficients();
    }

    /**
     * Private constructor.
     */
    private NumericalUtility() {
        //not called
    }
}
