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
import io.github.malapert.jwcs.proj.exception.JWcsError;
import io.github.malapert.jwcs.proj.exception.ProjectionException;
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.HALF_PI;
import static io.github.malapert.jwcs.utility.NumericalUtility.TWO_PI;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.util.FastMath;

/**
 * Conversion of intermediate world coordinates (or projection plane
 * coordinates) to celestial coordinates (\u03B1, \u03B4).
 *
 * <p>AbstractProjection plane is given and must be computed from pixels 
 * coordinates by applying a linear transformation.
 *
 * <p>The conversion is organized by a pipeline doing the following steps:
 * <ul>
 * <li>computes the projection plane coordinates (x,y) to native spherical coordinates
 * (\u03D5, \u03B8) through a spherical projection</li>
 * <li>computes the native spherical coordinates (\u03D5, \u03B8) to celestial spherical
 * coordinates through a spherical coordinate rotation.</li>
 * </ul>
 *
 * <p>Ref : "Representations of celestial coordinates in FITS", Calabretta, M.R.,
 * and Greisen, E.W., (2002), Astronomy and Astrophysics, 395, 1077-1122.
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public abstract class AbstractProjection {

    /**
     * Default native latitude of the celestial pole (\u03B8<sub>p</sub>) sets to {@link NumericalUtility#HALF_PI}.
     */
    public static final double DEFAULT_THETAP = HALF_PI;

    /**
     * Default native longitude of the celestial pole (\u03D5<sub>p</sub>) sets to 0.     
     */
    public static final double DEFAULT_PHIP = 0;
    /**
     * Logger.
     */
    private static final Logger LOG = Logger.getLogger(AbstractProjection.class.getName());

    /**
     * Native longitude in radians of the celestial pole for \u03B4<sub>0</sub>
     * &lt; \u03B8<sub>0</sub>.
     */
    protected static final double LONPOLE_PI = FastMath.PI;

    /**
     * Native longitude in radians of the celestial pole for \u03B4<sub>0</sub>
     * &ge; \u03B8<sub>0</sub>.
     */
    protected static final double LONPOLE_0 = 0;

    /**
     * Celestial longitude \u03B1<sub>0</sub> in radians of the fiducial point.
     */
    private double crval1;
    /**
     * Celestial longitude \u03B4<sub>0</sub> in radians of the fiducial point.
     */
    private double crval2;
    /**
     * Initializes the native longitude of the celestial pole (\u03D5<sub>p</sub>) to {@link AbstractProjection#DEFAULT_PHIP}. 
     */
    private double phip = DEFAULT_PHIP;
    /**
     * Initializes the native latitude of the celestial pole (\u03B8<sub>p</sub>) to {@link NumericalUtility#HALF_PI}. 
     */
    private double thetap = HALF_PI;    
    /**
     * Celestial longitude and latitude of the native pole (\u03B1<sub>p</sub>, \u03B4<sub>p</sub>).
     */
    private double[] coordNativePole;

    /**
     * Creates an instance of projection by given sky position coordinates.
     * 
     * <p>Initializes {@link AbstractProjection#crval1} with crval1, {@link AbstractProjection#crval2} with crval2 and {@link AbstractProjection#thetap} to {@link AbstractProjection#DEFAULT_THETAP}
     *
     * @param crval1 Celestial longitude (\u03B1<sub>0</sub>) in degrees of the ﬁducial point
     * @param crval2 Celestial latitude (\u03B4<sub>0</sub>) in degrees of the ﬁducial point
     */
    protected AbstractProjection(final double crval1, final double crval2) {
        this.crval1 = FastMath.toRadians(crval1);
        this.crval2 = FastMath.toRadians(crval2);
        LOG.log(Level.FINER, "INPUTS: crval1[deg]={0} crval2[deg]={1}", new Object[]{crval1,crval2});        
        LOG.log(Level.FINEST, "Theta_p[deg]={0}", new Object[]{FastMath.toDegrees(DEFAULT_THETAP)});
        setThetap(DEFAULT_THETAP);
    }

    /**
     * Computes the native spherical coordinates (\u03D5, \u03B8) from the projection plane
     * coordinates (x, y).
     *
     * @param x projection plane coordinate along X
     * @param y projection plane coordinate along Y
     * @return the native spherical coordinates (\u03D5, \u03B8) in radians
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException when
     * an error happens while the projection
     */
    protected abstract double[] project(final double x, final double y) throws ProjectionException;

    /**
     * Computes the projection plane coordinates (x, y) from the native spherical
     * coordinates (\u03D5, \u03B8).
     *
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param theta the native spherical coordinate (\u03B8) in radians along latitude
     * @return the projection plane coordinates
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException when
     * an error happens while the projection
     */
    protected abstract double[] projectInverse(final double phi, final double theta) throws ProjectionException;

    /**
     * Returns the native longitude of the fiducial point (\u03D5<sub>0</sub>) in radians.
     *
     * @return \u03D5<sub>0</sub> in radians.
     */
    public abstract double getPhi0();

    /**
     * Sets the native longitude of the ﬁducial point (\u03D5<sub>0</sub>) in radians.
     *
     * @param phio the native longitude in radians of the ﬁducial point (\u03D5<sub>0</sub>)
     */
    public abstract void setPhi0(final double phio);

    /**
     * Returns the native latitude of the ﬁducial point (\u03B8<sub>0</sub>) in radians.
     *
     * @return \u03B8<sub>0</sub> in radians.
     */
    public abstract double getTheta0();

    /**
     * Sets the native latitude in radians of the ﬁducial point (\u03B8<sub>0</sub>).
     *
     * @param theta0 the native latitude in radians of the ﬁducial point (\u03B8<sub>0</sub>)
     */
    public abstract void setTheta0(final double theta0);

    /**
     * Computes the default value for \u03D5<sub>p</sub>. 
     * 
     * <p>The default value of \u03D5<sub>p</sub> will be 
     * {@link AbstractProjection#LONPOLE_0} for \u03B4<sub>0</sub> &ge;
     * \u03B8<sub>0</sub> or {@link AbstractProjection#LONPOLE_PI} for
     * \u03B4<sub>0</sub> &lt; \u03B8<sub>0</sub>.
     *
     * @return \u03D5<sub>p</sub> in radians.
     */
    protected final double computeDefaultValueForPhip() {
        final double phi_p;
        if (NumericalUtility.equalValues(getCrval2(), getTheta0())) {
            phi_p = LONPOLE_0;
        } else if (getCrval2() > getTheta0()) {
            phi_p = LONPOLE_0;
        } else {
            phi_p = LONPOLE_PI;
        }
        return phi_p;
    }

    /**
     * Computes the celestial spherical coordinates (\u03B1, \u03B4) from the
     * native spherical coordinates (\u03D5, \u03B8). 
     * 
     * <p>The computation is performed by applying the spherical coordinate rotation.
     *
     * <p>general case:<br>
     * ------------<br>
     * \u03B1 = \u03B1<sub>p</sub> + arg(sin\u03B8cos\u03B4<sub>p</sub> -
     * cos\u03B8sin\u03B4<sub>p</sub>cos(\u03D5-\u03D5<sub>p</sub>),-cos\u03B8sin(\u03D5-\u03D5<sub>p</sub>))<br>
     * \u03B4 = asin(sin\u03B8sin\u03B4<sub>p</sub> +
     * cos\u03B8cos\u03B4<sub>p</sub>cos(\u03D5-\u03D5<sub>p</sub>))<br>
     * 
     * <p>Special cases:<br>
     * -------------<br>
     * if \u03B4<sub>p</sub> = HALF_PI<br>
     * \u03B1 = \u03B1<sub>p</sub> + \u03D5 - \u03D5<sub>p</sub> - PI<br>
     * \u03B4 = \u03B8<br>
     * if \u03B4<sub>p</sub> = -HALF_PI<br>
     * \u03B1 = \u03B1<sub>p</sub> - \u03D5 + \u03D5<sub>p</sub><br>
     * \u03B4 = -\u03B8<br>
     *     
     * @param phi Native longitude (\u03D5) in radians
     * @param theta Native latitude (\u03B8) in radians
     * @return Returns the celestial spherical coordinates (\u03B1, \u03B4) in
     * degrees
     */
    protected double[] computeCelestialSpherical(final double phi, final double theta) {
        double ra;
        double dec;
        LOG.log(Level.FINER, "INPUTS[deg]: (phi,theta)=({0},{1})", new Object[]{FastMath.toDegrees(phi),FastMath.toDegrees(theta)});        
        
        final double alphap = getCoordNativePole()[0];
        final double deltap = getCoordNativePole()[1];
        LOG.log(Level.FINEST, "CoordinateNativePole[deg]: (alphap,deltap)=({0},{1})", new Object[]{FastMath.toDegrees(alphap),FastMath.toDegrees(deltap)});        
        
        if (NumericalUtility.equalValues(deltap, HALF_PI)) {
            ra = alphap + phi - getPhip() - FastMath.PI;
            dec = theta;
        } else if (NumericalUtility.equalValues(deltap, -HALF_PI)) {
            ra = alphap - phi + getPhip();
            dec = -theta;
        } else {            
            ra = alphap + NumericalUtility.aatan2(-FastMath.cos(theta) * FastMath.sin(phi - phip),
                    FastMath.sin(theta) * FastMath.cos(deltap)
                    - FastMath.cos(theta) * FastMath.sin(deltap)
                    * FastMath.cos(phi - getPhip()),0);            

            dec = NumericalUtility.aasin(FastMath.sin(theta) * FastMath.sin(deltap)
                    + FastMath.cos(theta) * FastMath.cos(deltap)
                    * FastMath.cos(phi - getPhip()));
        }
        // convert ra, dec to degrees
        ra = FastMath.toDegrees(ra);
        dec = FastMath.toDegrees(dec);
        if (ra < 0) {
            ra += 360;
        }
        
        final double[] pos = {ra, dec};
        LOG.log(Level.FINER, "OUTPUTS[deg] pos=({0},{1})", new Object[]{FastMath.toDegrees(pos[0]),FastMath.toDegrees(pos[1])});
        return pos;
    }

    /**
     * Computes the native spherical coordinates (\u03D5, \u03B8) from 
     * the celestial spherical coordinates (\u03B1, \u03B4)  .
     * 
     * <p>The computation is performed by applying the inverse of the spherical 
     * coordinate rotation.<br>
     * \u03D5 = \u03D5<sub>p</sub> + arg(sin\u03B4cos\u03B4<sub>p</sub>-cos\u03B4sin\u03B4<sub>p</sub>cos(\u03B1-\u03B1<sub>p</sub>),-cos\u03B4sin(\u03B1-\u03B1<sub>p</sub>))<br>
     * \u03B8 = asin(sin\u03B4sin\u03B4<sub>p</sub>+cos\u03B4cos\u03B4<sub>p</sub>cos(\u03B1-\u03B1<sub>p</sub>))
     *
     * @param ra Celestial longitude in radians
     * @param dec Celestial latitude in radians
     * @return Returns native longitude and latitude in radians
     */
    protected double[] computeNativeSpherical(final double ra, final double dec) {
        LOG.log(Level.FINER, "INPUTS[deg]: (ra,dec)=({0},{1})", new Object[]{FastMath.toDegrees(ra),FastMath.toDegrees(dec)});        

        final double ra_p = getCoordNativePole()[0];
        final double dec_p = getCoordNativePole()[1];
        LOG.log(Level.FINEST, "CoordinateNativePole[deg]: (alphap,deltap)=({0},{1})", new Object[]{FastMath.toDegrees(ra_p),FastMath.toDegrees(dec_p)});                
        
        final double phi;
        final double theta;
        if (NumericalUtility.equalValues(dec_p, HALF_PI)) {
            phi = getPhip() + ra - ra_p;
            theta = dec;
        } else if (NumericalUtility.equalValues(dec_p, -HALF_PI)) {
            phi = getPhip() - ra + ra_p;
            theta = -dec;            
        } else {
            phi = getPhip() + NumericalUtility.aatan2(-FastMath.cos(dec) * FastMath.sin(ra - ra_p),
                    FastMath.sin(dec) * FastMath.cos(dec_p)
                    - FastMath.cos(dec) * FastMath.sin(dec_p)
                    * FastMath.cos(ra - ra_p), 0);
            theta = NumericalUtility.aasin(FastMath.sin(dec) * FastMath.sin(dec_p)
                     + FastMath.cos(dec) * FastMath.cos(dec_p)
                     * FastMath.cos(ra - ra_p));      
        }                   
        final double[] pos = {phi, theta};
        LOG.log(Level.FINER, "OUTPUTS[deg] (phi,theta)=({0},{1})", new Object[]{FastMath.toDegrees(phi),FastMath.toDegrees(theta)});
        return pos;
    }

    /**
     * Computes the celestial coordinates (\u03B1<sub>p</sub>,
     * \u03B4<sub>p</sub>) of the native pole
     * (\u03B4<sub>p</sub>=\u03B8<sub>p</sub>).
     *
     * <p>Projections such as the cylindricals and conics for which
     * (\u03D5<sub>0</sub>, \u03B8<sub>0</sub>) = (0, HALF_PI) are handled by
     * providing formulae to compute (\u03B1<sub>p</sub>, \u03B4<sub>p</sub>)
     * from (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     *
     * @param phi_p Native longitude (\u03D5<sub>p</sub>) in radians of the
     * celestial pole
     * @return Celestial (\u03B1<sub>p</sub>, \u03B4<sub>p</sub>) longitude and
     * latitude in radians of the native pole
     */
    protected double[] computeCoordNativePole(final double phi_p) {

        if (NumericalUtility.equalValues(getPhi0(), 0)
                && NumericalUtility.equalValues(getTheta0(), HALF_PI)) {
            LOG.log(Level.FINEST,"No need to compute the coordinates of the native pole");
            return new double[]{getCrval1(), getCrval2()};
        }
        // native longitude of the celestial pole in radians
        final double deltap = computeLatitudeNativePole(phi_p);
        final double alphap = computeLongitudeNativePole(deltap, phi_p);
        return new double[]{alphap, deltap};               
    }
        
    /**
     * Computes the latitude of the celestial coordinates (\u03B1<sub>p</sub>) 
     * of the native pole (\u03B4<sub>p</sub>=\u03B8<sub>p</sub>).
     * 
     * @param phi_p Native longitude (\u03D5<sub>p</sub>) in radians of the
     * celestial pole
     * @return Celestial (\u03B1<sub>p</sub>) longitude in radians of the native
     * pole
     * @throws JWcsError No valid solution for thetap
     */
    private double computeLatitudeNativePole(final double phi_p) {
        final double deltap;
        if (NumericalUtility.equalValues(getTheta0(), 0.0d) && NumericalUtility.equalValues(getCrval2(), 0) && NumericalUtility.equalValues(FastMath.abs(phi_p - getPhi0()), HALF_PI)) {
            deltap = getThetap();
        } else {
            final double deltap_arg = NumericalUtility.aatan2(FastMath.sin(getTheta0()), FastMath.cos(getTheta0()) * FastMath.cos(phi_p - getPhi0()));
            final double deltap_acos = NumericalUtility.aacos(FastMath.sin(getCrval2()) / FastMath.sqrt(1 - FastMath.pow(FastMath.cos(getTheta0()), 2) * FastMath.pow(FastMath.sin(phi_p - getPhi0()), 2)));
            final double deltap1 = deltap_arg + deltap_acos;
            final double deltap2 = deltap_arg - deltap_acos;

            if (NumericalUtility.equalValues(getTheta0(), 0)
                    && NumericalUtility.equalValues(getCrval2(), 0)
                    && NumericalUtility.equalValues(FastMath.abs(getPhip() - getPhi0()), HALF_PI)) {
                deltap = getThetap();
            } else {
                deltap = findTheValidDeltap(deltap1, deltap2);
            }
        }
        return deltap;
    }
    
    /**
     * Finds the valid celestial latitude of the native pole \u1D5F<sub>p</sub>.
     * 
     * <p>References:<br>
     * -----------<br>
     * Representations of celestial coordinates in FITS, 
     * Calabretta. M.R., and Greisen, E.w., (2002) Astronomy and Astrophysics,
     * 395, 1077-1122. http://www.atnf.csiro.au/people/mcalabre/WCS/ccs.pdf     
     * 
     * <p>Valid solutions are ones that lie in range -90&deg; to +90&deg;
     * @param deltap1 first solution
     * @param deltap2 second solution
     * @return <code>deltap1</code> or <code>deltap2</code>
     * @throws JWcsError No valid solution for thetap
     */
    private double findTheValidDeltap(final double deltap1, final double deltap2) {
        final double deltap;
        final boolean isDeltap1InInterval = NumericalUtility.isInInterval(deltap1, -HALF_PI, HALF_PI);
        final boolean isDeltap2InInterval = NumericalUtility.isInInterval(deltap2, -HALF_PI, HALF_PI);
        if (isDeltap1InInterval && isDeltap2InInterval) {
            final double diff1 = FastMath.abs(deltap1 - getThetap());
            final double diff2 = FastMath.abs(deltap2 - getThetap());
            deltap = diff1 < diff2 ? deltap1 : deltap2;
        } else if (isDeltap1InInterval) {
            deltap = deltap1;
        } else if (isDeltap2InInterval) {
            deltap = deltap2;
        } else {
            throw new JWcsError("No valid solution for thetap");
        }        
        return deltap;
    }
    
    /**
     * Computes the longitude of the celestial coordinates (\u03B4<sub>p</sub>) 
     * of the native pole (\u03B4<sub>p</sub>=\u03B8<sub>p</sub>).
     *
     * @param deltap Celestial latitude (\u03B4<sub>p</sub>) in radians of the native pole
     * @param phi_p Native longitude (\u03D5<sub>p</sub>) in radians of the
     * celestial pole
     * @return Celestial (\u03B1<sub>p</sub>) longitude in radians of the native
     * pole
     */    
    private double computeLongitudeNativePole(final double deltap, final double phi_p) {
        final double alphap;
        if (NumericalUtility.equalValues(FastMath.abs(getCrval2()), HALF_PI)) {
            alphap = getCrval1();
        } else if (NumericalUtility.equalValues(deltap, HALF_PI)) {
            alphap = getCrval1() + phi_p - getPhi0() - FastMath.PI;
        } else if (NumericalUtility.equalValues(deltap, -HALF_PI)) {
            alphap = getCrval1() - phi_p + getPhi0();      
        } else {
            final double das = FastMath.sin(phi_p - getPhi0()) * FastMath.cos(getTheta0()) / FastMath.cos(getCrval2());
            final double dac = (FastMath.sin(getTheta0()) - FastMath.sin(deltap) * FastMath.sin(getCrval2())) / (FastMath.cos(deltap) * FastMath.cos(getCrval2()));
            alphap = getCrval1() - NumericalUtility.aatan2(das, dac);
        }        
        return alphap;
    }

    /**
     * Sets the native longitude in radians of the celestial pole (\u03D5<sub>p</sub>).
     *
     * @param phip The native longitude in radians of the celestial pole (\u03D5<sub>p</sub>)
     */
    public final void setPhip(final double phip) {
        this.phip = phip;
        this.coordNativePole = computeCoordNativePole(phip);
    }

    /**
     * Returns the native longitude in radians of the celestial pole (\u03D5<sub>p</sub>).
     *
     * @return the native longitude in radians of the celestial pole (\u03D5<sub>p</sub>)
     */
    public double getPhip() {
        return this.phip;
    }

    /**
     * Sets the native latitude in radians of the celestial pole (\u03B8<sub>p</sub>).
     *
     * @param thetap the native latitude in radians of the celestial pole (\u03B8<sub>p</sub>)
     */
    public final void setThetap(final double thetap) {
        this.thetap = thetap;
    }

    /**
     * Returns the native latitude in radians of the celestial pole (\u03B8<sub>p</sub>).
     *
     * @return the native latitude in radians of the celestial pole (\u03B8<sub>p</sub>)
     */
    public double getThetap() {
        return this.thetap;
    }

    /**
     * Returns phi between [-PI, PI].
     *
     * @param phi phi in radians
     * @return phi between [-PI, PI]
     */
    protected final double phiRange(final double phi) {
        double phiCorrect = phi % TWO_PI;
        if (phiCorrect > FastMath.PI) {
            phiCorrect -= TWO_PI;
        } else if (phiCorrect < -FastMath.PI) {
            phiCorrect += TWO_PI;
        }
        return phiCorrect;
    }

    /**
     * Computes the celestial spherical coordinates (\u03B1, \u03B4) from the projection plane
     * coordinates (\u03B4, \u03B8).
     *
     * @param x projection plane longitude (\u03D5) in radians
     * @param y projection plane latitude (\u03B8) in radians
     * @return the celestial spherical coordinates in degrees (\u03B1, \u03B4)
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException when
     * an error happens while the projection
     */
    public double[] projectionPlane2wcs(final double x, final double y) throws ProjectionException {
        LOG.log(Level.FINER, "INPUTS[Deg] (x,y)=({0},{1})", new Object[]{x,y});                                                                                                                                                
        final double[] pos = project(x, y);
        LOG.log(Level.FINER, "OUTPUTS[Deg] (phi,theta)=({0},{1})", new Object[]{FastMath.toDegrees(pos[0]),FastMath.toDegrees(pos[1])}); 
        return computeCelestialSpherical(pos[0], pos[1]);
    }

    /**
     * Computes the projection plane coordinates (\u03B4, \u03B8) from the celestial spherical
     * coordinates (\u03B1, \u03B4).
     *
     * @param ra right ascension in radians (\u03B1)
     * @param dec declination in radians (\u03B4)
     * @return the projection plane coordinates (\u03D5, \u03B8)
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException when
     * an error happens while the projection
     */
    public double[] wcs2projectionPlane(final double ra, final double dec) throws ProjectionException {
        final double raFixed = NumericalUtility.normalizeLongitude(ra);
        double[] nativeSpherical = computeNativeSpherical(raFixed, dec);
        nativeSpherical[0] = phiRange(nativeSpherical[0]);
        LOG.log(Level.FINER, "INPUTS[Deg] (phi,theta)=({0},{1})", new Object[]{FastMath.toDegrees(nativeSpherical[0]), FastMath.toDegrees(nativeSpherical[1])});
        final double[] coord = projectInverse(nativeSpherical[0], nativeSpherical[1]);
        LOG.log(Level.FINER, "OUTPUTS[Deg] (x,y)=({0},{1})", coord); 
        return coord;
    }

    /**
     * Returns the celestial longitude in radians of the ﬁducial point (\u03B1<sub>0</sub>).
     *
     * @return \u03B1<sub>0</sub> in radians
     */
    public double getCrval1() {
        return this.crval1;
    }

    /**
     * Returns the celestial latitude in radians of the ﬁducial point (\u03B4<sub>0</sub>). 
     *
     * @return  \u03B4<sub>0</sub> in radians
     */
    public double getCrval2() {
        return this.crval2;
    }
    
    /**
     * Sets the crval1 in radians.
     * @param crval1 the crval1 to set
     */
    public void setCrval1(final double crval1) {
        this.crval1 = crval1;
    }

    /**
     * Sets the crval2 in radians.
     * @param crval2 the crval2 to set
     */
    public void setCrval2(final double crval2) {
        this.crval2 = crval2;
    }    

    /**
     * Returns true if the given lat/lon point is visible in this projection.
     *
     * @param lon longitude in radians (\u03B1).
     * @param lat latitude in radians (\u03B4).
     * @return True when the point is visible otherwise False.
     */
    public abstract boolean inside(final double lon, final double lat);
    
    /**
     * Checks if the line is visible.
     * @param pos1 First point of the line
     * @param pos2 Last point of the line
     * @return Returns True when the line is visible otherwise False
     */
    public abstract boolean isLineToDraw(final double[] pos1, final double[] pos2); 

    /**
     * Returns the projection's name.
     *
     * @return the projection's name
     */
    public abstract String getName();

    /**
     * Returns the projection's family.
     *
     * @return the projection's family
     */
    public abstract String getNameFamily();

    /**
     * Returns the projection's description.
     *
     * @return the projection's description
     */
    public abstract String getDescription();

    /**
     * Returns the celestial longitude and latitude of the native 
     * pole (\u03B1<sub>p</sub>, \u03B4<sub>p</sub>).
     * 
     * @return the coordNativePole (\u03B1<sub>p</sub>, \u03B4<sub>p</sub>)
     */
    protected double[] getCoordNativePole() {
        return coordNativePole.clone();
    }
    
    /**
     * Returns the logger of the projection family.
     * 
     * @return the logger
     */
    public abstract Logger getLogger();

    /**
     * Returns the projection parameters for a specific projection.
     * 
     * @return the projection parameters.
     */
    public abstract ProjectionParameter[] getProjectionParameters();
    
    /**
     * The ProjectionParameter object deFINERs few metadata about a projection parameter.
     * 
     * <p>This object is used in the GUI to display the projection parameter. It defines :
     * <ul>
     * <li>The name of the parameter.
     * <li>The PV keyword related to the name.
     * <li>The valid interval of the parameter
     * <li>The default value of the parameter to display in the GUI.    
     * </ul>
     * 
     * <p>To define an undefined value, the value is set to Double.POSITIVE_INFINITY 
     * for positive number and Double.NEGATIVE_INFINITY for negative number.
     */
    public static class ProjectionParameter {
        
        /**
         * Name of the parameter.
         */
        private final String name;
        /**
         * PV keyword related to the name of the parameter.
         */
        private final String pvName;
        /**
         * valid interval of the parameter.
         */
        private double[] validInterval;
        /**
         * Default value of the parameter.
         */
        private final double defaultValue;
        
        /**
         * Creates a new parameter.
         * 
         * @param name its name
         * @param pvName its related PV keyword in FITS
         * @param validInterval its valid interval
         * @param defaultValue  its default value.
         */
        public ProjectionParameter(final String name, final String pvName, final double[] validInterval, final double defaultValue) {
            this.name = name;
            this.pvName = pvName;
            setValidInterval(validInterval);
            this.defaultValue = defaultValue;
        }
        
        /**
         * Returns the name of the parameter.
         * 
         * @return the name of the parameter
         */
        public String getName() {
            return this.name;
        }
        
        /**
         * Returns the related keyword to PV.
         * 
         * @return the related keyword to PV
         */
        public String getPvName() {
            return this.pvName;
        }
        
        /**
         * Returns the valid interval of the parameter.
         * 
         * @return the valid interval of the parameter
         */
        public double[] getValidInterval() {
            return this.validInterval.clone();
        }
        
        /**
         * Returns the default value of the parameter.
         * 
         * @return the default value
         */
        public double getDefaultValue() {
            return this.defaultValue;
        } 

        /**
         * Sets the validity interval.
         * @param validInterval the validInterval to set
         */
        protected final void setValidInterval(final double[] validInterval) {
            if (validInterval == null) {
                this.validInterval = new double[0];
            } else {
                this.validInterval = Arrays.copyOf(validInterval, validInterval.length);
            }
        }
        
    }    
}
