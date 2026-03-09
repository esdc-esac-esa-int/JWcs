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
 * Cylindrical perspective.
 * 
 * <p>The sphere is projected onto a cylinder of radius <code>\u03BB</code> 
 * spherical radii from points in the equatorial plane of the native system at a
 * distance <code>\u03BC</code> spherical radii measured from the center of the 
 * sphere in the direction opposite the projected surface.
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class CYP extends AbstractCylindricalProjection {
    
    /**
     * Projection's name.
     */
    private static final String NAME_PROJECTION = "Cylindrical perspective";
    
    /**
     * Projection's description.
     */
    private static final String DESCRIPTION = "\u03BC=%s \u03BB=%s";

    /**
     * Default value for \u03BC.
     */
    public static final double DEFAULT_MU = 1;

    /**
     * Default value for \u03BB.
     */
    public static final double DEFAULT_LAMBDA = FastMath.sqrt(2)/2;
    
    /**
     * \u03BC: distance in spherical radii from the center of the sphere to the equatorial plane of the native system. 
     */
    private double mu;
    /**
     * \u03BB: radius of the cylinder in spherical radii.
     */
    private double lambda;

   /**
     * Constructs a CYP projection based on the default celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     *
     * <p>\u03BC is set to {@link CYP#DEFAULT_MU}.
     * \u03BB is set to {@link CYP#DEFAULT_LAMBDA}.
     * 
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException Lambda must be &gt; 0 or Mu must not be -lambda
     */    
    public CYP() throws BadProjectionParameterException {
        this(FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_THETA0));
    }
    
   /**
     * Constructs a CYP projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     *
     * <p>\u03BC is set to {@link CYP#DEFAULT_MU}.
     * \u03BB is set to {@link CYP#DEFAULT_LAMBDA}.
     * 
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException Lambda must be &gt; 0 or Mu must not be -lambda
     */
    public CYP(final double crval1, final double crval2) throws BadProjectionParameterException {
        this(crval1, crval2, DEFAULT_MU, DEFAULT_LAMBDA);
    }
    
    /**
     * Constructs a CYP projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) and \u03BC and \u03BB. 
     *     
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @param mu \u03BC distance measured from the center of the sphere in the direction opposite the projected surface
     * @param lambda \u03BB radius of the cylinder
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException Lambda must be &gt; 0 or Mu must not be -lambda
     * @see <a href="http://www.atnf.csiro.au/people/mcalabre/WCS/ccs.pdf">Representations of celestial coordinates in FITS, chapter 5.2.1</a>
     */
    public CYP(final double crval1, final double crval2, final double mu, final double lambda) throws BadProjectionParameterException {
        super(crval1, crval2);
        this.mu = mu;
        this.lambda = lambda;
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2,mu,lambda)=({0},{1},{2},{3})", new Object[]{crval1, crval2, mu, lambda});        
        checkParameters(mu, lambda);
    }
    
    /**
     * Checks projection parameters.
     * 
     * <p>Sets \u03BB = 1 when \u03BB &lt; 0 or \u03BC = -\u03BB.
     * @param mu value to check
     * @param lambda value to check
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException Lambda must be &gt; 0 or Mu must not be -lambda
     */
    protected final void checkParameters(final double mu, final double lambda) throws BadProjectionParameterException {
        if (getLambda() < 0 || NumericalUtility.equalValues(getLambda(), 0)) {
            throw new BadProjectionParameterException(this, "Lambda must be > 0");
        }
        if (NumericalUtility.equalValues(getMu(), -getLambda())) {
            throw new BadProjectionParameterException(this, "Mu must not be -lambda");
        }              
    }    

    @Override
    public double[] project(final double x, final double y) {
        final double xr = FastMath.toRadians(x);
        final double yr = FastMath.toRadians(y);
        final double phi = xr / getLambda();        
        final double eta = yr / (getMu() + getLambda());
        final double theta = NumericalUtility.aatan2(eta, 1) + NumericalUtility.aasin(getMu() * eta / FastMath.sqrt(FastMath.pow(eta, 2) + 1));
        return new double[] { phi, theta};
    }

    @Override
    public double[] projectInverse(final double phi, final double theta) throws PixelBeyondProjectionException {
        final double x = getLambda() * phi;
        final double ctheta = FastMath.cos(theta);
        if(NumericalUtility.equalValues(getMu(), -ctheta)) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(phi), FastMath.toDegrees(theta), false);
        }
        final double y = (getMu()+getLambda())/(getMu() + ctheta) * FastMath.sin(theta);
        return new double[] { FastMath.toDegrees(x), FastMath.toDegrees(y)};
    }

    /**
     * Returns \u03BC in spherical radii.
     * 
     * @return the mu
     */
    public final double getMu() {
        return mu;
    }

    /**
     * Returns \u03BB.
     * 
     * @return the lambda
     */
    public final double getLambda() {
        return lambda;
    }

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return String.format(DESCRIPTION, NumericalUtility.round(this.getMu()), NumericalUtility.round(this.getLambda()));
    }
    
    @Override
    public ProjectionParameter[] getProjectionParameters() {
        final ProjectionParameter p1 = new ProjectionParameter("\u03BC", AbstractJWcs.PV21, new double[]{Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}, 0);
        final ProjectionParameter p2 = new ProjectionParameter("\u03BB", AbstractJWcs.PV22, new double[]{0, Double.POSITIVE_INFINITY}, 1);
        return new ProjectionParameter[]{p1,p2};        
    }        

    /**
     * Sets \u03BC in spherical radii.
     * @param mu the mu to set
     */
    private void setMu(final double mu) {
        this.mu = mu;
    }

    /**
     * Sets \u03BB in spherical radii.
     * @param lambda the lambda to set
     */
    private void setLambda(final double lambda) {
        this.lambda = lambda;
    }
    
    /**
     * Sets the projection parameters.
     * @param mu mu
     * @param lambda lambda
     * @throws BadProjectionParameterException Lambda must be &gt; 0 or Mu must not be -lambda
     */
    public void setProjectionParameters(final double mu, final double lambda) throws BadProjectionParameterException {
        checkParameters(mu, lambda);
        setMu(mu);
        setLambda(lambda);
    }

}
