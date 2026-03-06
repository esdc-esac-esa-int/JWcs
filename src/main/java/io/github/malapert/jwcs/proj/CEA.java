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

import io.github.malapert.jwcs.AbstractJWcs;
import io.github.malapert.jwcs.proj.exception.BadProjectionParameterException;
import io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException;
import io.github.malapert.jwcs.utility.NumericalUtility;
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * The cylindrical equal area projection.
 * 
 * <p>Reference: "Representations of celestial coordinates in FITS", 
 * M. R. Calabretta and E. W. Greisen
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class CEA extends AbstractCylindricalProjection {

    /**
     * Projection's name.
     */
    private final static String NAME_PROJECTION = "Cylindrical equal area";

    /**
     * Projection's description.
     */
    private final static String DESCRIPTION = "\u03BB=%s";

    /**
     * \u03BB Scaling parameter.
     */
    private double lambda;

    /**
     * Default value for \u03BB.
     */
    private final static int DEFAULT_VALUE = 1;
    
    /**
     * Constructs a CEA based on the default celestial longitude and latitude of the
     * fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     * @throws BadProjectionParameterException When projection parameters are wrong
     */    
    public CEA() throws BadProjectionParameterException {
        this(FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_THETA0));
    }

    /**
     * Constructs a CEA based on the celestial longitude and latitude of the
     * fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     *
     * <p>\u03BB is set to {@link CEA#DEFAULT_VALUE}.
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @throws BadProjectionParameterException When projection parameters are wrong
     */
    public CEA(final double crval1, final double crval2) throws BadProjectionParameterException {
        this(crval1, crval2, DEFAULT_VALUE);
    }

    /**
     * Constructs a CEA based on the celestial longitude and latitude of the
     * fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) and \u03BB.
     * 
     * <p>\u03BB is set by the FITS keyword PV<code>nbAxis</code>_1.
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @param lambda \u03BB dimensionless.
     * @throws BadProjectionParameterException lambda not in ]0,1]
     */
    public CEA(final double crval1, final double crval2, final double lambda) throws BadProjectionParameterException {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2,lambda)=({0},{1},{2})", new Object[]{crval1,crval2,lambda});                        
        setLambda(lambda);
    }
    
    /**
     * Checks lambda parameter.
     * @param lambda \u03BB dimensionless
     * @throws BadProjectionParameterException lambda not in ]0,1]
     */
    private void checkParameter(final double lambda) throws BadProjectionParameterException {
        if (NumericalUtility.equal(lambda, 0) || lambda < 0 || lambda > 1.0) {
            throw new BadProjectionParameterException(this,"lambda =" + lambda + " - lambda outside the range ]0,1]");
        }        
    }

    /**
     * Computes the native spherical coordinates (\u03D5, \u03B8) from the projection plane
     * coordinates (x, y).
     * 
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes \u03D5 = x</li>
     * <li>computes \u03B8 = asin(\u03BB * y)</li>
     * </ul>
     * 
     * @param x projection plane coordinate along X
     * @param y projection plane coordinate along Y
     * @return the native spherical coordinates (\u03D5, \u03B8) in radians
     * @throws io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException No valid solution for (x,y)
     */     
    @Override
    protected double[] project(final double x, final double y) throws PixelBeyondProjectionException {
        final double xr = FastMath.toRadians(x);
        final double yr = FastMath.toRadians(y);
        final double phi = xr;
        final double arg = getLambda() * yr;
        final double theta = NumericalUtility.aasin(arg);
        if(Double.isNaN(theta)) {
            throw new PixelBeyondProjectionException(this, x, y, true);
        }
        return new double[] { phi, theta};
    }

    /**
     * Computes the projection plane coordinates (x, y) from the native spherical
     * coordinates (\u03D5, \u03B8).
     *
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes x = \u03D5</li>
     * <li>computes y = sin\u03B8 / \u03BB</li>
     * </ul>
     * 
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param theta the native spherical coordinate (\u03B8) in radians along latitude
     * @return the projection plane coordinates
     */    
    @Override
    protected double[] projectInverse(final double phi, final double theta) {
        final double x = FastMath.toDegrees(phi);
        final double y = FastMath.toDegrees(FastMath.sin(theta) / getLambda());
        return new double[] { x, y};
    }

    /**
     * Returns \u03BB, the scaling parameter.
     *
     * @return the lambda
     */
    public double getLambda() {
        return lambda;
    }

    /**
     * Sets \u03BB, the scaling parameter.
     *
     * @param lambda the lambda to set
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException lambda not in ]0,1]
     */
    public final void setLambda(final double lambda) throws BadProjectionParameterException {
        checkParameter(lambda);
        this.lambda = lambda;
    }

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return String.format(DESCRIPTION, NumericalUtility.round(this.lambda));
    }
    
    @Override
    public ProjectionParameter[] getProjectionParameters() {
        final ProjectionParameter p1 = new ProjectionParameter("\u03BB", AbstractJWcs.PV21, new double[]{0, 1}, 1);
        return new ProjectionParameter[]{p1};        
    }    

}
