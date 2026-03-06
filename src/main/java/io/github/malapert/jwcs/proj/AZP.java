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
import io.github.malapert.jwcs.proj.exception.MathematicalSolutionException;
import io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException;
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.HALF_PI;
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * Zenithal perspective.
 *
 * <p>Zenithal (azimuthal) perspective projections are generated from a point and
 * carried through the sphere to the plane of projection.
 *
 * @see <a href="http://www.atnf.csiro.au/people/mcalabre/WCS/ccs.pdf">
 * "Representations of celestial coordinates in FITS", M. R. Calabretta and
 * E. W. Greisen - page 10</a>
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public final class AZP extends AbstractZenithalProjection {

    /**
     * Projection's name.
     */
    private final static String NAME_PROJECTION = "Zenithal perspective";
    
    /**
     * Projection's description.
     */
    private final static String DESCRIPTION = "\u03BC=%s \u0263=%s";

    /**
     * \u0263 is the angle between the camera's optical axis and the line to the center of the planet.
     */
    private double gamma;
    /**
     * \u03BC is the distance from the center of the sphere to the source of projection.
     * \u03BC increases in the direction away from the plane of projection.
     */
    private double mu;
    /**
     * Default value for \u0263 and \u03BC.
     */
    private final static double DEFAULT_VALUE = 0;

    /**
     * Creates a new AZC projection based on the default celestial longitude and
     * latitude of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     * 
     * <p>\u03BC and \u0263 are set to {@link AZP#DEFAULT_VALUE}.
     *
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException Gamma must be different +/- 90°
     */    
    public AZP() throws BadProjectionParameterException {
        this(FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_THETA0));
    }
    
    /**
     * Creates a new AZC projection based on the celestial longitude and
     * latitude of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     * 
     * <p>\u03BC and \u0263 are set to {@link AZP#DEFAULT_VALUE}.
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the fiducial point
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException Gamma must be different +/- 90°
     */
    public AZP(final double crval1, final double crval2) throws BadProjectionParameterException {
        this(crval1, crval2, DEFAULT_VALUE, DEFAULT_VALUE);
    }

    /**
     * Creates a new AZP projection based on the celestial longitude and
     * latitude of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) and \u03BC and \u0263.
     * 
     * <p><img alt="AZP projection" src="doc-files/azpProjection.png">
     * <br>The picture shows alternate geometries of slant zenithal perspective 
     * projections with \u0263=2 and \u03BC=30&deg;
     * 
     * <p>\u03BC is in degrees while \u0263 is set in radii.
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the fiducial point
     * @param gamma \u03BC Angle in degrees between the camera's optical axis and the line to the center of the planet
     * @param mu \u0263 Distance in radii from the center of the sphere to the source of projectionPV<code>nbAxis</code>_1 in radians
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException Gamma must be different +/- 90°
     */
    public AZP(final double crval1, final double crval2, final double mu, final double gamma) throws BadProjectionParameterException {
        super(crval1, crval2);
        setGamma(FastMath.toRadians(gamma));
        setMu(mu);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2,mu,gamma)=({0},{1},{2},{3})", new Object[]{crval1,crval2,mu,gamma});        
    }
    
    /**
     * Checks \u0263 and \u03BC parameters.
     * 
     * @param mu value to check
     * @param gamma value to check
     * @throws BadProjectionParameterException \u03BC must be != +/- HALF_PI and \u0263 != -1 
     */
    private void checkParameters(final double mu, final double gamma) throws BadProjectionParameterException {       
        checkParameterMu(mu);
        checkParameterGamma(gamma);
    }
    
    /**
     * Checks the validity of \u03BC.
     * @param mu value to check
     * @throws BadProjectionParameterException \u03BC must be != -1
     */
    private void checkParameterMu(final double mu) throws BadProjectionParameterException {
        if (NumericalUtility.equal(mu, -1)) {
            throw new BadProjectionParameterException(this, "\u0263 must be != -1");
        }
    }
    
    /**
     * Checks the validity of \u0263.
     * @param gamma value to check
     * @throws BadProjectionParameterException \u0263 must be != +/-HALF_PI
     */
    private void checkParameterGamma(final double gamma) throws BadProjectionParameterException {
        if(NumericalUtility.equal(FastMath.abs(gamma), HALF_PI)) {
            throw new BadProjectionParameterException(this, "\u03BC must be != +/-HALF_PI");
        }        
    }
    
    /**
     * Computes the plane coordinate along X based on parameter projection.
     * 
     * <p>X coordinate is computed as <code>xr = x</code>
     * 
     * @param x plane coordinate along X in radians.
     * @return the plane coordinate along X with application of parameter projection
     */
    private double computeXr(final double x) {
        return x;
    }
    
    /**
     * Computes the plane coordinate along Y based on parameter projection.
     * 
     * <p>Y coordinate is computed as <code>yr = y * cos\u03BC</code>
     * 
     * @param y plane coordinate along Y in radians.
     * @return the plane coordinate along Y with application of parameter projection
     */    
    private double computeYr(final double y) {
        return y * FastMath.cos(getGamma());
    }       
    
    /**
     * Computes \u03C1.
     * 
     * <p>\u03C1 is computed as:
     * <br><code>\u03C1 = radius / (\u03BC + 1 + y * tan\u0263)</code>
     * 
     * @param y plane coordinate in radians along Y
     * @param radius Radius
     * @return rho
     * @throws BadProjectionParameterException getMu() + 1 + y * FastMath.tan(getGamma()) cannot be 0
     */
    private double computeRho(final double y, final double radius) throws BadProjectionParameterException {
        final double denom = getMu() + 1 + y * FastMath.tan(getGamma());
        if (NumericalUtility.equal(denom,0)) {
            throw new BadProjectionParameterException(this,"(mu,gamma) = (" + getMu() + ", " + getGamma()+"). (mu + 1) + y * tan(gamma) must be !=0");
        }    
        return radius / denom;        
    }    
    
    /**
     * Computes theta in radians.
     * 
     * <p>The algorithm to compute \u03B8 is the following:
     * <ul>
     * <li>computes \u03C1 : {@link AZP#computeRho(double, double) }</li>
     * <li>computes \u03C9 : asin(\u03BC * \u03C1 / sqrt(\u03C1<sup>2</sup> + 1)</li>
     * <li>computes \u0471 : arg(\u03C1, 1)</li>
     * <li>selects the right \u03B8 : {@link AZP#findTheValidSolution(double, double) }
     * for |\u03C9| &lt; 1 else {@link AZP#findTheSolutionNearestNorthPole(double, double) }</li>
     * </ul>     
     * 
     * @param x plane coordinate in radians along X
     * @param y plane coordinate in radians along Y
     * @param radius radius
     * @return the \u03B8 native spherical coordinates in radians
     * @throws BadProjectionParameterException When projection parameters are wrong
     * @throws PixelBeyondProjectionException When (x,y) has no solution
     */
    private double computeTheta(final double x, final double y, final double radius) throws BadProjectionParameterException, PixelBeyondProjectionException {   
        final double rho = computeRho(y, radius);
        
        //computes omega
        final double val = getMu()*rho/FastMath.sqrt(FastMath.pow(rho, 2)+1);
        final double omega = NumericalUtility.aasin(val);
        if (Double.isNaN(omega)) {
            throw new PixelBeyondProjectionException(this, x, y, true);
        }
        
        //computes psi
        final double psi = NumericalUtility.aatan2(1.0, rho);
        if(Double.isNaN(psi)) {
            throw new PixelBeyondProjectionException(this, x, y, true);            
        }
        
        //computes the valid theta solution
        final double theta1 = NumericalUtility.normalizeLatitude(psi - omega);
        final double theta2 = NumericalUtility.normalizeLatitude(psi + omega + FastMath.PI);        
        final double theta;        
        if(FastMath.abs(getMu()) < 1) {
            try {
                theta = findTheValidSolution(theta1, theta2);
            } catch (MathematicalSolutionException ex) {
                throw new PixelBeyondProjectionException(this, x, y, ex.getMessage(), true);
            }
        } else {
            theta = findTheSolutionNearestNorthPole(theta1, theta2);
        }
        return theta;
    }
        
    /**
     * Finds the right solution among theta1 and theta2.
     * 
     * <p>The right solution is this one which is the nearest from the north pole.
     * 
     * @param theta1 First solution
     * @param theta2 Second solution
     * @return theta1 or theta2
     */
    private double findTheSolutionNearestNorthPole(final double theta1, final double theta2) {
        final double theta;
        final double diffTheta1With90 = FastMath.abs(theta1-HALF_PI);
        final double diffTheta2With90 = FastMath.abs(theta2-HALF_PI);
        if (diffTheta1With90 > diffTheta2With90) {
            theta = theta2; 
        } else {
            theta = theta1;
        }    
        return theta;
    }
    
    /**
     * Finds the valid solution among the two solutions.
     * 
     * <p>Only one theta is valid. A valid theta is a theta in which the value
     * is included in [-HalfPi, HalfPi].
     * 
     * @param theta1 First solution
     * @param theta2 Second solution
     * @return theta1 or theta2
     * @throws MathematicalSolutionException No valid solution
     */
    private double findTheValidSolution(final double theta1, final double theta2) {
        final double theta;
        final boolean isTheta1Valid = NumericalUtility.isInInterval(theta1, -HALF_PI, HALF_PI);
        final boolean isTheta2Valid = NumericalUtility.isInInterval(theta2, -HALF_PI, HALF_PI);
        if(isTheta1Valid && isTheta2Valid) {
            throw new MathematicalSolutionException("No valid solution");
        } else if (isTheta1Valid) {
            theta = theta1;
        } else if (isTheta2Valid) {
            theta = theta2;
        } else {
            throw new MathematicalSolutionException("No valid solution");
        }
        return theta;
    }

    /**
     * Computes the native spherical coordinates (\u03D5, \u03B8) from the projection plane
     * coordinates (x, y).
     * 
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes radius : {@link AbstractZenithalProjection#computeRadius(double, double) }</li>
     * <li>computes \u03D5 : {@link AbstractZenithalProjection#computePhi(double, double, double) }</li>
     * <li>computes \u03B8 : {@link AZP#computeTheta(double, double, double) }</li>
     * </ul>
     *
     * @param x projection plane coordinate along X
     * @param y projection plane coordinate along Y
     * @return the native spherical coordinates (\u03D5, \u03B8) in radians
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException getMu() + 1 + y * FastMath.tan(getGamma()) cannot be 0
     * @throws io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException No valid solution for (x,y)
     */     
    @Override
    public double[] project(final double x, final double y) throws BadProjectionParameterException, PixelBeyondProjectionException {
        final double xr = computeXr(FastMath.toRadians(x));
        final double yr = computeYr(FastMath.toRadians(y));
        final double r = computeRadius(xr, yr);
        final double theta = computeTheta(xr, yr, r);
        final double phi = computePhi(xr, yr, r);
        return new double[] { phi, theta};
    }

    /**
     * Computes the projection plane coordinates (x, y) from the native spherical
     * coordinates (\u03D5, \u03B8).
     *
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes radius : {@link AZP#computeRadiusFrom(double, double) }</li>
     * <li>computes x : {@link AbstractZenithalProjection#computeX(double, double) }</li>
     * <li>computes y : {@link AbstractZenithalProjection#computeY(double, double) }</li>
     * </ul>
     * 
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param theta the native spherical coordinate (\u03B8) in radians along latitude
     * @return the projection plane coordinates
     * @throws io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException \u03B8 is beyond the projection
     */     
    @Override
    public double[] projectInverse(final double phi, final double theta) throws PixelBeyondProjectionException {
        final double r = computeRadiusFrom(phi, theta);
        final double x = computeX(r, phi);
        final double y = computeY(r / FastMath.cos(getGamma()), phi);
        return new double[] { FastMath.toDegrees(x), FastMath.toDegrees(y)};
    }
    
    /**
     * Computes the radius from (\u03D5, \u03B8).
     * 
     * <p>To compute (x,y), first we need to compute R, which is given by :
     * <br>R = (\u03BC + 1) * cos\u03B8 / (\u03BC + sin\u03B8 + cos\u03B8 * cos\u03D5 * tan\u0263)
     * 
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param theta the native spherical coordinate (\u03B8) in radians along latitude
     * @return the radius
     * @throws PixelBeyondProjectionException \u03B8 is beyond the projection
     */
    private double computeRadiusFrom(final double phi, final double theta) throws PixelBeyondProjectionException {
        final double denom = getMu() + FastMath.sin(theta) + FastMath.cos(theta) * FastMath.cos(phi) * FastMath.tan(getGamma());        
        if (!isVisible(theta, denom)) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(phi), FastMath.toDegrees(theta), false);
        }

        return (getMu() + 1) * FastMath.cos(theta) / denom;          
    }
    
    /**
     * Computes if theta is beyond the limb.
     * @param theta theta
     * @param denom denom
     * @return false when theta is beyond the limb
     */
    private boolean isVisible(final double theta, final double denom) {
        final double thetax;
        if (NumericalUtility.equal(getMu(), 0)) {
            thetax = 0;
        } else if (FastMath.abs(getMu()) > 1) {
            thetax = NumericalUtility.aasin(-1.0d / getMu());
        } else {
            thetax = NumericalUtility.aasin(-getMu());
        }        
        final boolean result;
        result = !(NumericalUtility.equal(denom, 0) || theta < thetax);        
        return result;
    }
    
    @Override
    public boolean inside(final double lon, final double lat) {
        final double raFixed = NumericalUtility.normalizeLongitude(lon);
        double[] nativeSpherical = computeNativeSpherical(raFixed, lat);
        nativeSpherical[0] = phiRange(nativeSpherical[0]);
        final double denom = getMu() + FastMath.sin(nativeSpherical[1]) + FastMath.cos(nativeSpherical[1]) * FastMath.cos(nativeSpherical[0]) * FastMath.tan(getGamma());                
        return isVisible(nativeSpherical[1], denom);
    }    

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }
    
    @Override
    public String getDescription() {
        return String.format(DESCRIPTION, this.getMu(), NumericalUtility.round(FastMath.toDegrees(this.getGamma())));
    }

    @Override
    public ProjectionParameter[] getProjectionParameters() {
        final ProjectionParameter p1 = new ProjectionParameter("\u03BC", AbstractJWcs.PV21, new double[]{Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}, 0);
        final ProjectionParameter p2 = new ProjectionParameter("\u0263", AbstractJWcs.PV22, new double[]{0, 360}, 0);
        return new ProjectionParameter[]{p1,p2};        
    }

    /**
     * Returns \u0263 in radians.
     * @return the gamma
     */
    public double getGamma() {      
        return gamma;
    }

    /**
     * Sets \u0263 in radians.
     * @param gamma the gamma to set
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException Gamma must be different +/- 90°
     */
    public void setGamma(final double gamma) throws BadProjectionParameterException {
        checkParameterGamma(gamma);
        this.gamma = gamma;
    }

    /**
     * Returns \u03BC.
     * @return the mu
     */
    public double getMu() {
        return mu;
    }

    /**
     * Sets \u03BC.
     * @param mu the mu to set
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException mu must be != -1
     */
    public void setMu(final double mu) throws BadProjectionParameterException {
        checkParameterMu(mu);
        this.mu = mu;
    }
    
    /**
     * Sets the projection parameters.
     * @param mu mu
     * @param gamma gamma
     * @throws BadProjectionParameterException \u03BC must be != +/- HALF_PI and \u0263 != -1 
     */
    public void setProjectionParameters(final double mu, final double gamma) throws BadProjectionParameterException {
        checkParameters(mu, gamma);
        setMu(mu);
        setGamma(gamma);
    }
}
