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
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.HALF_PI;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.util.FastMath;

/**
 * In conic projections the sphere is thought to be projected onto
 * the surface of a cone which is then opened out. 
 * 
 * <p>The native coordinate system is chosen so that the poles are coincident with
 * the axis of the cone. Native meridians are then projected as uniformly spaced
 * rays that intersect at a point (either directly or by extrapolation), and 
 * parallels are projected as equiangular arcs of concentric circles.
 * 
 * <p><img alt="View of the zenithal projection" src="doc-files/conicProjection.png">
 * 
 * <p>Two-standard conic projections are characterized by two
 * latitudes, theta1 and theta2, whose parallels are projected at their true
 * length. In the conic perspective projection these are the latitudes at which 
 * the cone intersects the sphere. 
 * <br>One-standard conic projections have 
 * theta1 = theta2 and the cone is tangent to the sphere as shown. 
 * Since conics are designed to minimize distortion in the regions between the
 * two standard parallels they are constructed so that the point on the prime 
 * meridian mid-way between the two standard parallels maps to the reference 
 * point
 * 
 * <p>Ref : "Representations of celestial coordinates in FITS", Calabretta, M.R., 
 * and Greisen, E.W., (2002), Astronomy and Astrophysics, 395, 1077-1122. - p11
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public abstract class AbstractConicProjection extends AbstractProjection {
    /**
     * Logger.
     */
    protected static final Logger LOG = Logger.getLogger(AbstractConicProjection.class.getName());
    /**
     * AbstractProjection name.
     */
    public static final String NAME = "Conic projections";
    /**
     * thetaA = (theta1 + theta2) / 2 in radians.
     */
    private final double thetaA;
    /**
     * eta = abs(theta1 - theta2) / 2 in radians.
     */
    private double eta;
    /**
     * Native longitude value in radians for conic projection.
     */
    protected static final double DEFAULT_PHI0 = 0;
    /**
     * Native longitude in radians of the ﬁducial point for the conic
     * projection.
     */
    private double phio;
    /**
     * Native latitude in radians of the ﬁducial point for the conic projection.
     */
    private double theta0;
    
    /**
     * thetaA - eta.
     */
    private double theta1;
    
    /**
     * thetaA + eta.
     */
    private double theta2;
    
    /**
     * Creates a new conic projection.
     * 
     * @param crval1 Celestial longitude in degrees of the ﬁducial point
     * @param crval2 Celestial latitude in degrees of the ﬁducial point
     * @param theta_a (theta1 + theta2) / 2 in degrees
     * @param eta abs(theta1 - theta2) / 2 in degrees
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException Each angle must be -90&le;theta1,theta2&le;90"
     */
    protected AbstractConicProjection(final double crval1, final double crval2, final double theta_a, final double eta) throws BadProjectionParameterException {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[deg] (crval1,crval2,theta_a,eta) = ({0},{1},{2},{3})", new Object[]{crval1, crval2, theta_a, eta});
        this.thetaA = FastMath.toRadians(theta_a);
        this.eta = FastMath.toRadians(eta);        
        this.theta1 = this.thetaA - this.eta;
        this.theta2 = this.thetaA + this.eta;
        LOG.log(Level.FINEST, "(theta1,theta2)[deg]=({0},{1})", new Object[]{FastMath.toDegrees(this.theta1),FastMath.toDegrees(this.theta2)}); 
        checkParameters(theta1, theta2);
        setPhi0(DEFAULT_PHI0);
        setTheta0(this.thetaA);
        setPhip(computeDefaultValueForPhip());
        LOG.log(Level.FINEST, "(phi0,theta0)[DEG]=({0},{1})", new Object[]{FastMath.toDegrees(DEFAULT_PHI0), FastMath.toDegrees(this.thetaA)});
        LOG.log(Level.FINEST, "phip[deg]={0}", FastMath.toDegrees(computeDefaultValueForPhip()));         
    }
    
    /**
     * Checks \u03B8<sub>1</sub>,\u03B8<sub>2</sub> parameters.
     * 
     * @param theta1 \u03B8<sub>1</sub> in radians
     * @param theta2 \u03B8<sub>2</sub> in radians
     * @throws BadProjectionParameterException When (theta1,theta2) not in range [-90,90]
     */
    private void checkParameters(final double theta1, final double theta2) throws BadProjectionParameterException {
        final boolean inRangeTheta1 = NumericalUtility.isInInterval(theta1, -HALF_PI, HALF_PI);
        final boolean inRangeTheta2 = NumericalUtility.isInInterval(theta2, -HALF_PI, HALF_PI);
        if(!inRangeTheta1 || !inRangeTheta2) {
            throw new BadProjectionParameterException(this,"(theta1,theta2). Each angle must be -90<=theta1,theta2<=90");
        }        
    }
    
    
    
    @Override
    public String getNameFamily() {
        return NAME;
    }         
    
    /**
     * Computes the native spherical coordinate (\u03D5) in radians along longitude.
     * 
     * @param x projection plane coordinate along X
     * @param y projection plane coordinate along Y
     * @param r_theta radius
     * @param y0 y0
     * @param c constant of the cone
     * @return native spherical coordinate (\u03D5) in radians along longitude
     */
    protected double computePhi(final double x, final double y, final double r_theta, final double y0, final double c) {
        return NumericalUtility.equal(r_theta, 0) ? 0 : NumericalUtility.aatan2(x/r_theta, (y0-y)/r_theta)/c;
    }
    
    /**
     * Computes the projection plane coordinate along X.
     * 
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param r_theta radius
     * @param c constant of the cone
     * @return the projection plane coordinate along X
     */
    protected double computeX(final double phi, final double r_theta, final double c) {
        return r_theta * FastMath.sin(c*phi);
    }

    /**
     * Computes the projection plane coordinate along Y.
     * 
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param r_theta radius
     * @param c constant of the cone
     * @param y0 y0 in radians
     * @return the projection plane coordinate along Y
     */
    protected double computeY(final double phi, final double r_theta, final double c, final double y0) {
        return -r_theta * FastMath.cos(c*phi) + y0;
    }    

    @Override
    public final double getPhi0() {
        return phio;
    }

    @Override
    public final double getTheta0() {
        return theta0;
    }  
    
    @Override
    public final void setPhi0(final double phio) {
        this.phio = phio;
    }

    @Override
    public final void setTheta0(final double theta0) {
        this.theta0 = theta0;
    }    

    /**
     * Returns thetaA in radians.
     * @return the thetaA
     */
    protected final double getThetaA() {
        return thetaA;
    }

    /**
     * Returns eta in radians.
     * 
     * @return the eta
     */
    protected final double getEta() {
        return eta;
    }
    
    /**
     * Sets eta in radians.
     * 
     * @param eta eta
     */
    protected void setEta(final double eta) {
        this.eta = eta;
    }      
    
    @Override
    public boolean inside(final double lon, final double lat) {     
        final double angle = NumericalUtility.distAngle(new double[]{getCrval1(), getCrval2()}, new double[]{lon, lat});
        LOG.log(Level.FINER, "(lont,lat,distAngle)[deg] = ({0},{1}) {2}", new Object[]{FastMath.toDegrees(lon), FastMath.toDegrees(lat), angle});
        return NumericalUtility.equal(angle, HALF_PI) || angle <= HALF_PI;
    }   
    
    @Override
    public boolean isLineToDraw(final double[] pos1, final double[] pos2) {
        LOG.log(Level.FINER, "True");
        return true;
    }    

    @Override
    public final Logger getLogger() {
        return LOG;
    }       

    /**
     * Returns theta1.
     * @return the theta1
     */
    protected final double getTheta1() {
        return theta1;
    }

    /**
     * Sets theta1.
     * @param theta1 the theta1 to set
     */
    protected final void setTheta1(final double theta1) {
        this.theta1 = theta1;
    }

    /**
     * Returns theta2.
     * @return the theta2
     */
    protected final double getTheta2() {
        return theta2;
    }

    /**
     * Sets theta2.
     * @param theta2 the theta2 to set
     */
    protected final void setTheta2(final double theta2) {
        this.theta2 = theta2;
    }
}
