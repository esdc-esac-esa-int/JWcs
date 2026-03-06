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

import io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException;
import io.github.malapert.jwcs.utility.GammaFunction;
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.HALF_PI;
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * Mollweide's.
 *
 * <p>In Mollweide's pseudocylindrical projection17, the meridians are projected as
 * ellipses that correctly divide the equator and the parallels are spaced so as
 * to make the projection equal area
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class MOL extends AbstractCylindricalProjection {

    /**
     * Projection's name.
     */
    private final static String NAME_PROJECTION = "Mollweide’s";

    /**
     * Projection's description.
     */
    private final static String DESCRIPTION = "no limits";

    /**
     * Default maximum iteration for the iterative solution.
     */
    public final static int DEFAULT_MAX_ITER = 1000;

    /**
     * Maximum iteration for the iterative solution.
     */
    private int maxIter;
    
    /**
     * GammaFunction function to solve.
     */
    private final GammaFunction gammaFunction;

    /**
     * Constructs a MOL projection based on the default celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     */    
    public MOL() {
        this(FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_THETA0));
    }
    
    /**
     * Constructs a MOL projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     */
    public MOL(final double crval1, final double crval2) {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2)=({0},{1})", new Object[]{crval1, crval2});
        setMaxIter(DEFAULT_MAX_ITER);
        this.gammaFunction = new GammaFunction();
    }

    @Override
    protected double[] project(final double x, final double y) throws PixelBeyondProjectionException {
        final double xr = FastMath.toRadians(x);
        final double yr = FastMath.toRadians(y);        
        final double[] phis = computePhiAndS(xr, yr);
        final double phi = phis[0];
        final double s = phis[1];
        final double theta = computeTheta(xr, yr, s);
        return new double[] { phi, theta};
    }

    /**
     * Computes the native spherical coordinate (\u03B8) in radians along
     * latitude.
     *
     * @param xr projection plane coordinate along X in radians
     * @param yr projection plane coordinate along Y in radians
     * @param s s
     * @return the native spherical coordinate (\u03B8) in radians along
     * latitude
     * @throws PixelBeyondProjectionException Solution not defined
     */
    private double computeTheta(final double xr, final double yr, final double s) throws PixelBeyondProjectionException {
        double z = yr / FastMath.sqrt(2);
        if (NumericalUtility.equal(FastMath.abs(z), 1)) {
            z = (z < 0.0 ? -1.0 : 1.0) + s * yr / FastMath.PI;
        } else if (FastMath.abs(1) > 1) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(xr), FastMath.toDegrees(yr), true);
        } else {
            z = NumericalUtility.aasin(z) / HALF_PI + s * yr / FastMath.PI;
        }
        if (NumericalUtility.equal(FastMath.abs(z), 1)) {
            z = z < 0.0 ? -1.0 : 1.0;
        } else if (FastMath.abs(1) > 1) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(xr), FastMath.toDegrees(yr), true);
        }
        final double theta = NumericalUtility.aasin(z);
        if (Double.isNaN(theta)) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(xr), FastMath.toDegrees(yr), true);
        }
        return theta;
    }

    /**
     * Computes the native spherical coordinate (\u03D5) in radians along
     * longitude and s.
     *
     * @param xr projection plane coordinate along X in radians
     * @param yr projection plane coordinate along Y in radians     
     * @return the native spherical coordinate (\u03D5) in radians along
     * longitude and s
     * @throws PixelBeyondProjectionException Solution not defined
     */
    private double[] computePhiAndS(final double xr, final double yr) throws PixelBeyondProjectionException {
        final double tol = 1.0e-12;        
        double s = 2 - FastMath.pow(yr, 2);
        final double phi;
        if (s <= tol) {
            if (s < -tol) {
                throw new PixelBeyondProjectionException(this, FastMath.toDegrees(xr), FastMath.toDegrees(yr), true);
            }
            s = 0.0;
            if (FastMath.abs(xr) > tol) {
                throw new PixelBeyondProjectionException(this, FastMath.toDegrees(xr), FastMath.toDegrees(yr), true);
            }
            phi = 0;
        } else {
            s = FastMath.sqrt(s);
            phi = HALF_PI * xr / s;
        }
        return new double[]{phi, s};
    }

    @Override
    protected double[] projectInverse(final double phi, final double theta) {
        final double gamma = computeGamma(theta);
        final double x = FastMath.toDegrees((FastMath.sqrt(2.0d) / HALF_PI) * phi * FastMath.cos(gamma));
        final double y = FastMath.toDegrees(FastMath.sqrt(2.0d) * FastMath.sin(gamma));
        return new double[] { x, y};
    }

    /**
     * Computes gamma by an iterative solution. 
     * 
     * <p>Solves <code>v - PI*sin(theta) + sin(v) = 0</code><br>
     * with <code>gamma = 0.5 * v</code>
     *
     * @param theta the native spherical coordinate (\u03B8) in radians along
     * latitude
     * @return gamma
     * @see GammaFunction
     */
    private double computeGamma(final double theta) {
        this.gammaFunction.setTheta(theta);
        return NumericalUtility.computeFunctionSolution(this.getMaxIter(), this.gammaFunction, -FastMath.PI, FastMath.PI) * 0.5;
    }

    /**
     * Returns the maximal number of iterations of the approximative solution of
     * the inverse projection.
     *
     * @return the number maximum of iteration
     */
    public int getMaxIter() {
        return maxIter;
    }

    /**
     * Sets the maximal number of iterations.
     *
     * @param maxIter the maximum number of iteration to set
     */
    public final void setMaxIter(final int maxIter) {
        this.maxIter = maxIter;
    }

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return DESCRIPTION;
    }   
}
