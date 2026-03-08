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
 * Conic equal area.
 *
 * <p>The standard parallels in Alber's conic equal area projection are projected
 * as concentric arcs at their true length and separated so that the area
 * between them is the same as the corresponding area on the sphere. The other
 * parallels are then drawn as concentric arcs spaced so as to preserve the
 * area.
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class COE extends AbstractConicProjection {

    /**
     * Projection's name.
     */
    private static final String NAME_PROJECTION = "Conic equal area";
    
    /**
     * Projection's description.
     */
    private static final String DESCRIPTION = "\u03B8a=%s \u03B7=%s";
        
    /**
     * Constant in radians.
     *
     * <p>This constant is equal to sin\u03B8<sub>a</sub>     
     */
    private final double y0;
    
    /**
     * Constant in radians.
     * 
     * <p>This constant is equal to sin\u03B8<sub>1</sub>+sin\u03B8<sub>2</sub>
     */
    private final double gamma;
    
    /**
     * Constant of the cone in radians.
     * 
     * <p>This constant is equal to 0.5*\u0263
     */
    private final double c;  
    
    /**
     * Constructs a COE projection based on the default celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) and 
     * 03B8<sub>a</sub> = 45 and \u03B7 = 25.
     *
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException When projection parameters are wrong
     */    
    public COE() throws BadProjectionParameterException {
        this(FastMath.toDegrees(AbstractConicProjection.DEFAULT_PHI0), 45, 45, 25);        
    }
    
   /**
     * Constructs a COE projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) and \u03B8<sub>a</sub> and \u03B7.
     *
     * <p>\u03B8<sub>a</sub> is set by the FITS keyword PV<code>nbAxis</code>_1 in degrees.
     * \u03B7 is set by the FITS keyword PV<code>nbAxis</code>_2 in degrees.
     * 
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @param theta_a \u03B8<sub>a</sub> in degrees and defined as \u03B8<sub>a</sub>=(\u03B8<sub>1</sub>+\u03B8<sub>2</sub>)/2
     * @param eta \u03B7 in degrees and defined as \u03B7=|\u03B8<sub>1</sub>-\u03B8<sub>2</sub>|/2
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException Projection parameters: sin(theta1) + sin(theta2) = 0 are not allowed
     */
    public COE(final double crval1, final double crval2, final double theta_a, final double eta) throws BadProjectionParameterException {
        super(crval1, crval2, theta_a, eta);  
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2,theta_a,eta)=({0},{1},{2},{3})", new Object[]{crval1,crval2,theta_a,eta});        
        gamma = FastMath.sin(getTheta1()) + FastMath.sin(getTheta2());
        checkParameters(gamma);
        c = gamma * 0.5;       
        y0 = FastMath.sqrt(1.0d + FastMath.sin(getTheta1()) * FastMath.sin(getTheta2()) - gamma * FastMath.sin((getTheta1()+getTheta2())*0.5)) / c;
    }
    
    /**
     * Checks parameters.
     * @param gamma parameter to check
     * @throws BadProjectionParameterException Projection parameters: sin(theta1) + sin(theta2) = 0 are not allowed
     */
    private void checkParameters(final double gamma) throws BadProjectionParameterException {
        if(NumericalUtility.equal(gamma, 0)) {
            throw new BadProjectionParameterException(this, "Projection parameters: sin(theta1) + sin(theta2) = 0 are not allowed");
        }
    }

    /**
     * Computes the native spherical coordinates (\u03D5, \u03B8) from the projection plane
     * coordinates (x, y).
     * 
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes radius : sign(\u03B8<sub>a</sub>) * sqrt(x<sup>2</sup> + ({@link COE#y0} - y)<sup>2</sup>)</li>
     * <li>computes \u03D5 : {@link AbstractConicProjection#computePhi(double, double, double, double, double) }</li>
     * <li>computes \u03B8 : asin(1 / {@link COE#gamma} + sin\u03B8<sub>1</sub> * sin\u03B8<sub>2</sub> / {@link COE#gamma} - {@link COE#gamma} * (0.5 * radius)<sup>2</sup>)</li>
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
        final double r_theta = FastMath.signum(getThetaA()) * FastMath.sqrt(FastMath.pow(xr, 2) + FastMath.pow(getY0() - yr, 2));
        final double phi = computePhi(xr, yr, r_theta, getY0(), getC());                   
        final double w = 1.0d / getGamma() + FastMath.sin(getTheta1()) * FastMath.sin(getTheta2()) / getGamma() - getGamma() * FastMath.pow(r_theta * 0.5, 2);
        final double theta = NumericalUtility.aasin(w);
        if (Double.isNaN(theta)) {
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
     * <li>computes radius : sqrt( 1 + sin\u03B8<sub>1</sub> * sin\u03B8<sub>2</sub> - {@link COE#gamma} * sin\u03B8 / {@link COE#c}
     * <li>computes x : {@link AbstractConicProjection#computeX(double, double, double) }</li>
     * <li>computes y : {@link AbstractConicProjection#computeY(double, double, double, double) } </li>
     * </ul>
     * 
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param theta the native spherical coordinate (\u03B8) in radians along latitude
     * @return the projection plane coordinates
     */  
    @Override
    protected double[] projectInverse(final double phi, final double theta) {
        LOG.log(Level.FINER, "INPUTS[Deg] (phi,theta)=({0},{1})", new Object[]{FastMath.toDegrees(phi),FastMath.toDegrees(theta)});                                                                
        final double r_theta = FastMath.sqrt(1.0d + FastMath.sin(getTheta1()) * FastMath.sin(getTheta2()) - getGamma() * FastMath.sin(theta)) / getC();      
        final double x = computeX(phi, r_theta, getC());
        final double y = computeY(phi, r_theta, getC(), getY0());
        final double[] coord = {FastMath.toDegrees(x), FastMath.toDegrees(y)};
        LOG.log(Level.FINER, "OUTPUTS[Deg] (x,y)=({0},{1})", coord);                                                        
        return coord;
    }
    
    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return String.format(DESCRIPTION, NumericalUtility.round(FastMath.toDegrees(this.getThetaA())), NumericalUtility.round(FastMath.toDegrees(this.getEta())));
    }   
    
    @Override
    public boolean inside(final double lon, final double lat) {
        LOG.log(Level.FINER, "True");        
        return true;
    }    
    
    @Override
    public ProjectionParameter[] getProjectionParameters() {
        final ProjectionParameter p1 = new ProjectionParameter("\u03B8a", AbstractJWcs.PV21, new double[]{-90, 90}, -45);
        final ProjectionParameter p2 = new ProjectionParameter("\u03B7", AbstractJWcs.PV22, new double[]{0, 90}, 0);
        return new ProjectionParameter[]{p1,p2};    
    }    

    /**
     * Returns gamma.
     * @return the gamma
     */
    private double getGamma() {
        return gamma;
    }

    /**
     * Returns the constant of the cone.
     * @return the c
     */
    private double getC() {
        return c;
    }

    /**
     * Returns y0.
     * @return the y0
     */
    private double getY0() {
        return y0;
    }

}
