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
import io.github.malapert.jwcs.proj.exception.JWcsError;
import io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException;
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.HALF_PI;
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * Slant zenithal perspective.
 *
 * <p>While the generalization of the AZP projection to tilted planes of projection
 * is useful for certain applications it does have a number of drawbacks, in
 * particular, unequal scaling at the reference point. 
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class SZP extends AbstractZenithalProjection {
    
    /**
     * Projection's name.
     */
    private final static String NAME_PROJECTION = "Slant zenithal perspective";
    
    /**
     * Projection's description.
     */
    private final static String DESCRIPTION = "\u03BC=%s \u03C6c=%s \u03B8c=%s";      

    /**
     * \u03BC : Distance in spherical radii from the center of the sphere to the source of the projection.
     */
    private double mu;
    /**
     * Intersection of the line PO with the sphere at the \u03D5<sub>c</sub> coordinate.
     */
    private double thetac;
    /**
     * Intersection of the line PO with the sphere at the \u03B8<sub>c</sub> coordinate.
     */    
    private double phic;
    
    /**
     * X coordinate of P.
     */
    private double xp;
    
    /**
     * Y coordinate of P.
     */
    private double yp;
    
    /**
     * Z coordinate of P.
     */
    private double zp;

    /**
     * Default value for \u03BC.
     */
    public final static double DEFAULT_VALUE_MU = 0;

    /**
     * Default value for \u03D5<sub>c</sub>.
     */
    public final static double DEFAULT_VALUE_PHIC = 0;

    /**
     * Default value for \u03B8<sub>c</sub>.
     */
    public final static double DEFAULT_VALUE_THETAC = 90;

   /**
     * Constructs a SZP projection based on the default celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     * 
     * <p>\u03D5<sub>c</sub> is set to {@link SZP#DEFAULT_VALUE_PHIC}.
     * \u03B8<sub>c</sub> is set to {@link SZP#DEFAULT_VALUE_THETAC}.
     * 
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException When projection parameters are wrong
     */    
    public SZP() throws BadProjectionParameterException {
        this(FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_THETA0));
    }
    
   /**
     * Constructs a SZP projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     * 
     * <p>\u03D5<sub>c</sub> is set to {@link SZP#DEFAULT_VALUE_PHIC}.
     * \u03B8<sub>c</sub> is set to {@link SZP#DEFAULT_VALUE_THETAC}.
     * 
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException When projection parameters are wrong
     */
    public SZP(final double crval1, final double crval2) throws BadProjectionParameterException {
        this(crval1, crval2, DEFAULT_VALUE_MU, DEFAULT_VALUE_THETAC, DEFAULT_VALUE_PHIC);
    }

   /**
     * Constructs a SZP projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     * 
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @param mu \u03BC parameter projection
     * @param phic \u03B8<sub>c</sub> parameter projection
     * @param thetac \u03D5<sub>c</sub> parameter projection
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException When projection parameters are wrong
     */
    public SZP(final double crval1, final double crval2, final double mu, final double phic, final double thetac) throws BadProjectionParameterException {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2,mu,phic,thetac)=({0},{1},{2},{3},{4})", new Object[]{crval1,crval2,mu,phic,thetac});                                                                                                                                                
        this.mu = mu;
        this.thetac = FastMath.toRadians(thetac);
        this.phic = FastMath.toRadians(phic);
        init();
        checkParameters(this.mu, this.phic, this.thetac);        
    }
    
    /**
     * Init xp,yp,zp.
     */
    private void init() {
        this.xp = -this.mu * FastMath.cos(this.thetac) * FastMath.sin(this.phic);
        this.yp = this.mu * FastMath.cos(this.thetac) * FastMath.cos(this.phic);
        this.zp = this.mu * FastMath.sin(this.thetac) + 1;          
    }

    /**
     * Check.
     * 
     * @param mu mu
     * @param phic phic
     * @param thetac thetac
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException When projection parameters are wrong
     * @throws JWcsError Non-standard phi0 or theta0 values
     */
    protected final void checkParameters(final double mu, final double phic, final double thetac) throws BadProjectionParameterException {
        if (!NumericalUtility.equal(getPhi0(), 0) || !NumericalUtility.equal(getTheta0(),HALF_PI)) {
            throw new JWcsError("Non-standard phi0 or theta0 values");
        }
        if (NumericalUtility.equal(this.zp, 0)) {
            throw new BadProjectionParameterException(this,"zp = 0. It must be !=0");
        }
    }

    /**
     * Computes the native spherical coordinates (\u03D5, \u03B8) from the projection plane
     * coordinates (x, y).
     * 
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes \u03B8 : {@link NumericalUtility#computeQuatraticSolution(double[]) }</li>
     * <li>computes \u03D5 : {@link AbstractZenithalProjection#computePhi(double, double, double) }</li>      
     * </ul>
     * 
     * @param x projection plane coordinate along X
     * @param y projection plane coordinate along Y
     * @return the native spherical coordinates (\u03D5, \u03B8) in radians     
     * @throws io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException No mathematical solution found   
     */      
    @Override
    public double[] project(final double x, final double y) throws PixelBeyondProjectionException {
        final double xr = FastMath.toRadians(x);
        final double yr = FastMath.toRadians(y);        
        final double X = xr;
        final double Y = yr;
        final double X1 = (X - xp) / zp;
        final double Y1 = (Y - yp) / zp;
        final double a = X1 * X1 + Y1 * Y1 + 1;
        final double b = (X1 * (X - X1) + Y1 * (Y - Y1))*2;
        final double c = (X - X1) * (X - X1) + (Y - Y1) * (Y - Y1) - 1;
        final double theta;
        try {
            theta = NumericalUtility.computeQuatraticSolution(new double[]{c,b,a});
        } catch (MathematicalSolutionException ex) {
            throw new PixelBeyondProjectionException(this, x, y, ex.getMessage(), true);
        }
        final double phi = computePhi(X - X1 * (1 - FastMath.sin(theta)), Y - Y1 * (1 - FastMath.sin(theta)), 1);
        return new double[] { phi, theta};
    }

    /**
     * Computes the projection plane coordinates (x, y) from the native spherical
     * coordinates (\u03D5, \u03B8).
     *
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes denom : zp - (1 - sin(\u03B8))</li>     
     * <li>computes x : (zp * cos(\u03B8) * sin(\u03D5) - x * (1 - sin(\u03B8)))/denom;</li>
     * <li>computes y : -(zp * cos(\u03B8) * cos(\u03D5) + yp * (1 - sin(\u03B8)))/denom</li>
     * </ul>
     * 
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param theta the native spherical coordinate (\u03B8) in radians along latitude
     * @return the projection plane coordinates
     * @throws io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException No valid solution for (\u03D5, \u03B8)
     */     
    @Override
    public double[] projectInverse(final double phi, final double theta) throws PixelBeyondProjectionException {
        final double denom = zp - (1 - FastMath.sin(theta));        
        if (!isVisible(phi, theta, denom)) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(phi), FastMath.toDegrees(theta), false);            
        }
        final double x = (zp * FastMath.cos(theta) * FastMath.sin(phi) - xp * (1 - FastMath.sin(theta)))/denom;
        final double y = -(zp * FastMath.cos(theta) * FastMath.cos(phi) + yp * (1 - FastMath.sin(theta)))/denom;
        return new double[] { FastMath.toDegrees(x), FastMath.toDegrees(y)};
    }  
    
    /**
     * Computes if theta is beyond the limb.
     * @param phi phi phi
     * @param theta theta 
     * @param denom denom
     * @return false when theta is beyond the limb
     */
    private boolean isVisible(final double phi, final double theta, final double denom) {
        if (!firstContstraintVisibility(phi, theta)) {
            return false;
        }
        if (!secondConstraintVisibility(theta)) {
            return false;
        }        
        return thirdConstraintVisibility(denom);
    }
    
    /**
     * Evaluates the constraint.
     * @param denom denom
     * @return true when it is visible
     */
    private boolean thirdConstraintVisibility(final double denom) {
        return !NumericalUtility.equal(denom, 0);     
    }
    
    /**
     * Evaluates the constraint.
     * @param theta theta
     * @return true when it is visible
     */    
    private boolean secondConstraintVisibility(final double theta) {
        final double thetaLimit = NumericalUtility.aasin(1 - this.zp);
        if (Double.isNaN(thetaLimit)) {
            return true;
        }
        return theta > thetaLimit;
    }
    
    /**
     * Evaluates the constraint.
     * @param phi phi
     * @param theta theta
     * @return true when it is visible
     */    
    private boolean firstContstraintVisibility(final double phi, final double theta) {
        final double rho = getMu() * FastMath.sin(getThetac());
        final double sigma = -getMu() * FastMath.cos(getThetac()) * FastMath.cos(phi - getPhic());
        final double omega = NumericalUtility.aasin(1.0d / FastMath.sqrt(FastMath.pow(rho, 2) + FastMath.pow(sigma, 2)));
        final double psi = NumericalUtility.aatan2(sigma, rho);
        double thetax1 = psi - omega;
        if (!Double.isNaN(thetax1)) {
            thetax1 = NumericalUtility.normalizeLatitude(thetax1)+1e-4;
        }
        double thetax2 = psi + omega + FastMath.PI;
        if (!Double.isNaN(thetax2)) {
            thetax2 = NumericalUtility.normalizeLatitude(thetax2)+1e-4;
        }        
        boolean result;
        if (NumericalUtility.isInInterval(thetax1, -HALF_PI, HALF_PI) && NumericalUtility.isInInterval(thetax2, -HALF_PI, HALF_PI)) {
            result = thetax1 < theta && theta < thetax2;
        } else if (NumericalUtility.isInInterval(thetax1, -HALF_PI, HALF_PI)) {
            result = theta > thetax1;
        } else if (NumericalUtility.isInInterval(thetax2, -HALF_PI, HALF_PI)) {
            result = theta > thetax2;
        } else {
            result = false;
        }
        return result;
    }
        

    @Override
    public boolean inside(final double lon, final double lat) {
        final double raFixed = NumericalUtility.normalizeLongitude(lon);
        final double[] nativeSpherical = computeNativeSpherical(raFixed, lat);
        nativeSpherical[0] = phiRange(nativeSpherical[0]);
        final double denom = zp - (1 - FastMath.sin(nativeSpherical[1])); 
        return isVisible(nativeSpherical[0], nativeSpherical[1], denom);
    }
    
    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return String.format(DESCRIPTION, NumericalUtility.round(this.getMu()), NumericalUtility.round(FastMath.toDegrees(this.getPhic())), NumericalUtility.round(FastMath.toDegrees(this.getThetac())));
    }
    
    @Override
    public ProjectionParameter[] getProjectionParameters() {
        final ProjectionParameter p1 = new ProjectionParameter("\u03BC", AbstractJWcs.PV21, new double[]{Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}, 0);
        final ProjectionParameter p2 = new ProjectionParameter("\u03C6c", AbstractJWcs.PV22, new double[]{0, 360}, 0);                
        final ProjectionParameter p3 = new ProjectionParameter("\u03B8c", AbstractJWcs.PV23, new double[]{0, 90}, 90);
        return new ProjectionParameter[]{p1,p2,p3};    
    }    

    /**
     * Returns mu.
     * @return the mu
     */
    public double getMu() {
        return mu;
    }

    /**
     * Returns thetac.
     * @return the thetac
     */
    public double getThetac() {
        return thetac;
    }

    /**
     * Returns phic.
     * @return the phic
     */
    public double getPhic() {
        return phic;
    }
    
    /**
     * Set projection parameters.
     * @param mu mu
     * @param phic phic
     * @param thetac thetac
     * @throws BadProjectionParameterException TODO
     */
    private void setProjectionParameters(final double mu, final double phic, final double thetac) throws BadProjectionParameterException {
        checkParameters(mu, phic, thetac);
        setMu(mu);
        setThetac(thetac);
        setPhic(phic);
    }

    /**
     * Sets mu.
     * @param mu the mu to set
     */
    private void setMu(final double mu) {
        this.mu = mu;
    }

    /**
     * Sets thetac in radians.
     * @param thetac the thetac to set
     */
    private void setThetac(final double thetac) {
        this.thetac = thetac;
    }

    /**
     * Sets phic in radians.
     * @param phic the phic to set
     */
    private void setPhic(final double phic) {
        this.phic = phic;
    }

}
