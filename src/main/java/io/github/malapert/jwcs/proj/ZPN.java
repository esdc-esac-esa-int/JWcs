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
package io.github.malapert.jwcs.proj;

import io.github.malapert.jwcs.proj.exception.BadProjectionParameterException;
import io.github.malapert.jwcs.proj.exception.MathematicalSolutionException;
import io.github.malapert.jwcs.proj.exception.JWcsError;
import io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException;
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.HALF_PI;
import java.util.Arrays;
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * Zenithal polynomial.
 *
 * <p>The zenithal polynomial projection, ZPN, generalizes the ARC projection by
 * adding polynomial terms up to a large degree in the zenith distance
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public final class ZPN extends AbstractZenithalProjection {

    /**
     * Projection's name.
     */
    private final static String NAME_PROJECTION = "Zenithal polynomial";

    /**
     * Projection's description.
     */
    private final static String DESCRIPTION = "poly=(%s)";

    /**
     * Default maximum iteration for iterative solution.
     */
    public final static int DEFAULT_MAX_ITER = 1000;

    /**
     * Maximum iteration for iterative solution.
     */
    private int maxIter;
    /**
     * Projection parameters.
     */
    private double[] pv;
    
    /**
     * Order of the polynomial function.
     */
    private final int n;

    /**
     * Creates a ZPN projection based on crval1, crval2 and the default projection
     * parameters and the pv = [0.050, 0.975, -0.807, 0.337, -0.065, 0.010, 0.003, -0.001].
     *
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException
     * when a parameter projection is wrong
     */    
    public ZPN() throws BadProjectionParameterException {
        this(FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_THETA0), new double[]{0.050, 0.975, -0.807, 0.337, -0.065, 0.010, 0.003, -0.001});
    }

    /**
     * Creates a ZPN projection based on crval1, crval2 and the projection
     * parameters.
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @param pv projection parameters
     * @throws
     * io.github.malapert.jwcs.proj.exception.BadProjectionParameterException
     * when a parameter projection is wrong
     */
    public ZPN(final double crval1, final double crval2, final double[] pv) throws BadProjectionParameterException {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2)=({0},{1} PV={2})", new Object[]{crval1, crval2, Arrays.toString(pv)});
        setPv(pv);
        setMaxIter(DEFAULT_MAX_ITER);
        check();
        this.n= NumericalUtility.getPolynomialOrder(pv);
    }


    /**
     * Checks validity of parameters.
     *
     * @throws JWcsError Non-standard PVi_1 and/or PVi_2 values
     * @throws ArrayIndexOutOfBoundsException Need at least 10 projection
     * parameters
     */
    private void check() {
        if (!NumericalUtility.equal(getPhi0(), 0) || !NumericalUtility.equal(getTheta0(), HALF_PI)) {
            throw new JWcsError("Non-standard PVi_1 and/or PVi_2 values");
        }
    }

    /**
     * Returns the number of maximum iteration for the iterative solution.
     *
     * @return the maximum number or iteration
     */
    public int getMaxIter() {
        return maxIter;
    }

    /**
     * Sets the number of maximum iteration for the iterative solution.
     *
     * @param maxIter the maxIter to set
     */
    public void setMaxIter(final int maxIter) {
        this.maxIter = maxIter;
    }

    /**
     * The highest pv coefficient not equal to 0.
     *
     * @return the n
     */
    private int getN() {
        return n;
    }
    
    /**
     * Computes a polynomial where orders are given by pv.
     *
     * @param x value
     * @param pv polynomial order
     * @return the value of the polynomial
     */
    private double polyEval(final double x, final double[] pv) {
        final int lastElt = pv.length - 1;
        double y = 0;
        final double result;
        for (int i = lastElt; i >= 0; i--) {
            y = y * x + pv[i];
        }
        result = y;

        return result;
    }

    /**
     * Computes the solution for a linear equation.
     *
     * @param f polynomial function
     * @return the solution for a linear equation
     */
    private double linearSolution(final Object f) {
        final double[] coefficients = NumericalUtility.getPolynomialCoefficients(f);
        return -coefficients[0] / coefficients[1];
    }

    /**
     * Solves the solution for a quadratic equation.
     *
     * @param f polynomial coefficients
     * @return the solution for a quadratic equation
     * @throws MathematicalSolutionException No mathematical solution found
     * @see NumericalUtility#computeQuatraticSolution
     */
    private double quadraticSolution(final Object f) throws MathematicalSolutionException {
        final double[] coeff = NumericalUtility.getPolynomialCoefficients(f);
        return NumericalUtility.computeQuatraticSolution(coeff);
    }

    /**
     * Solves a polynomial function.
     * @param f the polynomial function
     * @return the solution
     * @see NumericalUtility#computePolynomialSolution
     */
    private double polynomialSolution(final Object f) {
        final double[] coeff = NumericalUtility.getPolynomialCoefficients(f);
        final double result = coeff[0] > 0 ? 0 : NumericalUtility.computePolynomialSolution(this.getMaxIter(), f, 0, FastMath.PI);
        return result;
    }

    /**
     * Compute the solution for whatever polynomial equation.
     *
     * @param polynomialFunction polynomial function
     * @return the solution for whatever polynomial equation
     * @throws MathematicalSolutionException No mathematical solution found
     */
    private double computeSolution(final Object polynomialFunction) throws MathematicalSolutionException {
        final double result;
        switch (getN()) {
            case 1:
                result = linearSolution(polynomialFunction);
                break;
            case 2:
                result = quadraticSolution(polynomialFunction);
                break;
            default:
                result = polynomialSolution(polynomialFunction);
                break;
        }
        return result;
    }

    /**
     * Computes the native spherical coordinates (\u03D5, \u03B8) from the projection plane
     * coordinates (x, y).
     * 
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes radius : {@link ZPN#computeRadius(double, double) }</li>
     * <li>computes \u03D5 : {@link AbstractZenithalProjection#computePhi(double, double, double) }</li>      
     * <li>computes \u03B8 : HALF_PI - {@link ZPN#computeSolution(java.lang.Object) }</li>
     * </ul>
     * 
     * @param x projection plane coordinate along X
     * @param y projection plane coordinate along Y
     * @return the native spherical coordinates (\u03D5, \u03B8) in radians
     * @throws io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException No valid solution for (x,y)
     */     
    @Override
    protected double[] project(final double x, final double y) throws PixelBeyondProjectionException {
        try {
            final double xr = FastMath.toRadians(x);
            final double yr = FastMath.toRadians(y);
            final double r_theta = computeRadius(xr, yr);
            final double[] coeffPolynomial = getPv();
            coeffPolynomial[0] = coeffPolynomial[0] - r_theta;
            final Object polynomialFunction = NumericalUtility.createPolynomialFunction(coeffPolynomial);
            final double phi = computePhi(xr, yr, r_theta);
            final double theta = HALF_PI - computeSolution(polynomialFunction);
            final double[] pos = {phi, theta};
            return pos;
        } catch (MathematicalSolutionException ex) {
            throw new PixelBeyondProjectionException(this, x, y, ex.getMessage(), true);
        }
    }

    /**
     * Computes the projection plane coordinates (x, y) from the native spherical
     * coordinates (\u03D5, \u03B8).
     *
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes radius : {@link ZPN#polyEval(double, double[]) }</li>
     * <li>computes x : {@link AbstractZenithalProjection#computeX(double, double) }</li>
     * <li>computes y : {@link AbstractZenithalProjection#computeY(double, double) }</li>
     * </ul>
     * 
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param theta the native spherical coordinate (\u03B8) in radians along latitude
     * @return the projection plane coordinates
     */    
    @Override
    protected double[] projectInverse(final double phi, final double theta) {
        final double r_theta = FastMath.toDegrees(polyEval(HALF_PI - theta, getPv()));        
        final double x = computeX(r_theta, phi);
        final double y = computeY(r_theta, phi);
        final double[] coord = {x, y};
        return coord;
    }

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return String.format(DESCRIPTION, Arrays.toString(this.getPv()));
    }

    @Override
    public boolean inside(final double lon, final double lat) {
        return true;
    }

    @Override
    public ProjectionParameter[] getProjectionParameters() {
        return new ProjectionParameter[]{};
    }

    /**
     * Returns a copy of pv.
     *
     * @return the pv
     */
    public double[] getPv() {
        return pv.clone();
    }

    /**
     * Sets pv.
     *
     * @param pv the pv to set
     */
    public void setPv(final double[] pv) {
        if (pv == null) {
            this.pv = new double[0];
        } else {
            this.pv = Arrays.copyOf(pv, pv.length);
        }
    }

}
