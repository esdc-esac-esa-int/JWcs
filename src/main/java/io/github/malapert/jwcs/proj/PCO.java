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
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.HALF_PI;
import io.github.malapert.jwcs.utility.PcoFunction;
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * Polyconic.
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class PCO extends AbstractPolyConicProjection {

    /**
     * Projection's name.
     */
    private final static String NAME_PROJECTION = "Polyconic";

    /**
     * Projection's description.
     */
    private final static String DESCRIPTION = "no limits";

    /**
     * Default tolerance for the approximative solution of the inverse
     * projection.
     */
    public final static double DEFAULT_TOLERANCE = 1E-16;

    /**
     * Default maximum iterations for the approximative solution of the inverse
     * projection.
     */
    public final static int DEFAULT_MAX_ITER = 1000;

    /**
     * Number of iterations for the approximative solution.
     */
    private int maxIter;

    /**
     * Function to solve theta.
     */
    private final PcoFunction pcoFunction;

    /**
     * Constructs a PCO projection based on the default celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     */    
    public PCO() {
        this(FastMath.toDegrees(AbstractPolyConicProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractPolyConicProjection.DEFAULT_THETA0));
    }
    
    /**
     * Constructs a PCO projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     */
    public PCO(final double crval1, final double crval2) {
        super(crval1, crval2, 45);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2)=({0},{1},45)", new Object[]{crval1, crval2});
        setMaxIter(DEFAULT_MAX_ITER);
        this.pcoFunction = new PcoFunction();
    }


    /**
     * Returns the maximal number of iterations.
     *
     * @return the maxIter
     */
    public int getMaxIter() {
        return maxIter;
    }

    /**
     * Sets the maximal number of iterations.
     *
     * @param maxIter the maximum number of iteration
     */
    public final void setMaxIter(final int maxIter) {
        this.maxIter = maxIter;
    }

    @Override
    protected double[] project(final double x, final double y) throws PixelBeyondProjectionException {
        final double xr = FastMath.toRadians(x);
        final double yr = FastMath.toRadians(y);
        final double phi;
        final double theta;
        if (NumericalUtility.equal(yr, 0)) {
            phi = xr;
            theta = 0.0;
        } else if (NumericalUtility.equal(yr, HALF_PI)) {
            phi = 0.0;
            theta = yr < 0.0 ? -HALF_PI : HALF_PI;
        } else {
            final double[] position = computeIterativeSolution(xr, yr);
            phi = position[0];
            theta = position[1];
        }

        return new double[] { phi, theta};
    }

    @Override
    protected double[] projectInverse(final double phi, final double theta) {
        final double costhe = FastMath.cos(theta);
        final double sinthe = FastMath.sin(theta);
        final double a = phi * sinthe;
        double x;
        double y;
        if (NumericalUtility.equal(sinthe, 0.0)) {
            x = phi;
            y = 0.0;
        } else {
            final double cotthe = costhe / sinthe;
            x = cotthe * FastMath.sin(a);
            y = cotthe * (1.0 - FastMath.cos(a)) + theta;
        }
        return new double[] { FastMath.toDegrees(x), FastMath.toDegrees(y)};
    }
    
    /**
     * Computes the iterative solution of the {@link PcoFunction}  using 
     * bisection algorithm.
     * 
     * @param xr projection plane coordinate along X in radians
     * @param yr projection plane coordinate along Y in radians
     * @return an array representing in the order phi and theta
     */
    private double[] computeIterativeSolution(final double xr, final double yr) {
        final double min;
        final double max;        
        if (yr > 0.0) {
            min = 0;
            max = HALF_PI;
        } else {
            min = -HALF_PI;
            max = 0;
        }
        this.pcoFunction.set(xr, yr);
        final double theta = NumericalUtility.computeFunctionSolution(getMaxIter(), pcoFunction, min, max);  
        final double tanthe = FastMath.tan(theta);
        final double xp = 1 - (yr - theta) * tanthe;
        final double yp = xr * tanthe;
        final double phi;
        if (NumericalUtility.equal(xp, 0) && NumericalUtility.equal(yp, 0)) {
            phi = 0.0;
        } else {
            phi = NumericalUtility.aatan2(yp, xp) / FastMath.sin(theta);
        }
        return new double[]{phi, theta};
    }    

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return DESCRIPTION;
    }

    @Override
    public ProjectionParameter[] getProjectionParameters() {
        return new ProjectionParameter[]{};
    }

}
