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
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * Slant orthographic.
 *
 * <p>It represents the visual appearance of a sphere, e.g. a planet, when seen
 * from a great distance. 
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class SIN extends AbstractZenithalProjection {

    /**
     * Projection's name.
     */
    private final static String NAME_PROJECTION = "Slant orthographic";

    /**
     * Projection's description.
     */
    private final static String DESCRIPTION = "\u046F=%s \u03B7=%s";

    /**
     * Default value.
     */
    public final static double DEFAULT_VALUE = 0;

    /**
     * \u03BE is defined as \u03BE = cot\u03B8<sub>c</sub>sin\u03D5<sub>c</sub>.    
     */
    private double ksi;
    
    /**
     * \u03B7 is defined as \u03B7 = -cot\u03B8<sub>c</sub>cos\u03D5<sub>c</sub>.
     */
    private double eta;

   /**
     * Constructs a SIN projection based on the default celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) with default \u03BE,\u03B7 
     * 
     * <p>\u03BE,\u03B7 parameters are set to {@link SIN#DEFAULT_VALUE}.
     */    
    public SIN() {
        this(FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_THETA0));
    }
    
   /**
     * Constructs a SIN projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) with default \u03BE,\u03B7 
     * 
     * <p>\u03BE,\u03B7 parameters are set to {@link SIN#DEFAULT_VALUE}.
     * 
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     */
    public SIN(final double crval1, final double crval2) {
        this(crval1, crval2, DEFAULT_VALUE, DEFAULT_VALUE);
    }

   /**
     * Constructs a SIN projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) and \u03BE,\u03B7.      
     * 
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @param ksi \u03BE dimensionless
     * @param eta \u03B7 dimensionless
     */    
    public SIN(final double crval1, final double crval2, final double ksi, final double eta) {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2,ksi,eta)=({0},{1},{2},{3})", new Object[]{crval1,crval2,ksi,eta});                                                
        this.ksi = ksi;
        this.eta = eta;
    }

    @Override
    public double[] project(final double x, final double y) throws PixelBeyondProjectionException {
        final double xr = FastMath.toRadians(x);
        final double yr = FastMath.toRadians(y);
        final double phi;
        final double theta;
        if (NumericalUtility.equal(getKsi(), DEFAULT_VALUE) && NumericalUtility.equal(getEta(), DEFAULT_VALUE)) {
            final double r_theta = computeRadius(xr, yr);
            phi = computePhi(xr, yr, r_theta);
            
            //computes theta
            if(!NumericalUtility.isInInterval(r_theta,0 ,true ,1 ,false)) {
                throw new PixelBeyondProjectionException(this,x ,y , "r_theta must be < 1", true);
            }             
            theta = NumericalUtility.aacos(r_theta);
            
        } else {
            final double[] coeffFromReducedDiscrimant = computeCoeffOfReducedDiscrimant(xr, yr);
            try {
                theta = NumericalUtility.computeQuatraticSolution(coeffFromReducedDiscrimant);
            } catch (MathematicalSolutionException ex) {
                throw new PixelBeyondProjectionException(this, x, y, ex.getMessage(), true);
            }

            phi = NumericalUtility.aatan2(xr - getKsi() * (1 - FastMath.sin(theta)), -(yr - eta * (1 - FastMath.sin(theta))));
        }

        return new double[] { phi, theta};
    }
    
    /**
     * Computes if theta is beyond the limb.
     * @param phi phi phi
     * @param theta theta     
     * @return false when theta is beyond the limb
     */
    private boolean isVisible(final double phi, final double theta) {
        final double thetax = -FastMath.atan(ksi*FastMath.sin(phi)-eta*FastMath.cos(phi));
        return theta > thetax;
    }    
    
    /**
     * Computes the coefficients of the reduced discriminant.
     * @param xr projection plane coordinate along X in radians
     * @param yr projection plane coordinate along Y in radians
     * @return the coefficients of the reduced discrimant as (c,b,a)
     */
    private double[] computeCoeffOfReducedDiscrimant(final double xr, final double yr) {
        final double a = FastMath.pow(getKsi(), 2) + FastMath.pow(getEta(), 2) + 1;
        final double b = (getKsi() * (xr - getKsi()) + getEta() * (yr - getEta())) * 2;
        final double c = FastMath.pow(xr - getKsi(),2) + FastMath.pow(yr - getEta(),2) - 1;
        return new double[]{c,b,a};
    }

    @Override
    public double[] projectInverse(final double phi, final double theta) throws PixelBeyondProjectionException {
        if (!isVisible(phi, theta)) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(phi), FastMath.toDegrees(theta), false);
        }
        final double x = FastMath.cos(theta) * FastMath.sin(phi) + getKsi() * (1 - FastMath.sin(theta));
        final double y = -FastMath.cos(theta) * FastMath.cos(phi) + getEta() * (1 - FastMath.sin(theta));
        return new double[] { FastMath.toDegrees(x), FastMath.toDegrees(y)};
    }
    
    @Override
    public boolean inside(final double lon, final double lat) {
        final double raFixed = NumericalUtility.normalizeLongitude(lon);
        double[] nativeSpherical = computeNativeSpherical(raFixed, lat);
        nativeSpherical[0] = phiRange(nativeSpherical[0]);
        return isVisible(nativeSpherical[0], nativeSpherical[1]);
    }     

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return String.format(DESCRIPTION, NumericalUtility.round(this.getKsi()), NumericalUtility.round(this.getEta()));
    }
    
    @Override
    public ProjectionParameter[] getProjectionParameters() {
        final ProjectionParameter p1 = new ProjectionParameter("\u03BE", AbstractJWcs.PV21, new double[]{Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}, 0);
        final ProjectionParameter p2 = new ProjectionParameter("\u03B7", AbstractJWcs.PV22, new double[]{Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}, 0);
        return new ProjectionParameter[]{p1,p2};    
    }    

    /**
     * Returns \u03BE dimensionless.
     * @return the ksi
     */
    public double getKsi() {
        return ksi;
    }

    /**
     * Sets \u03BE dimensionless.
     * @param ksi the ksi to set
     */
    public void setKsi(final double ksi) {
        this.ksi = ksi;
    }

    /**
     * Returns \u03B7 dimensionless.
     * @return the eta
     */
    public double getEta() {
        return eta;
    }

    /**
     * Sets \u03B7 dimensionless.
     * @param eta the eta to set
     */
    public void setEta(final double eta) {
        this.eta = eta;
    }        
}
