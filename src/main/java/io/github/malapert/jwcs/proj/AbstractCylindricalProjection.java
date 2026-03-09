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

import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.util.FastMath;

/**
 * Cylindrical projections are so named because the surface of
 * projection is a cylinder. 
 * 
 * <p>The native coordinate system is chosen to have its polar axis coincident 
 * with the axis of the cylinder. Meridians and parallels are mapped onto a 
 * rectangular graticule.
 * 
 * <p><img alt="View of the zenithal projection" src="doc-files/cylindricalProjection.png">
 * <br>The sphere is projected onto a cylinder of radius \u03BB spherical radii 
 * from points in the equatorial plane of the native system at a distance \u03BC
 * spherical radii measured from the center of the sphere in the direction opposite
 * the projected surface
 * 
 * <p>Ref : "Representations of celestial coordinates in FITS", Calabretta, M.R., 
 * and Greisen, E.W., (2002), Astronomy and Astrophysics, 395, 1077-1122. - p15
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public abstract class AbstractCylindricalProjection extends AbstractProjection {
    /**
     * Logger.
     */
    protected static final Logger LOG = Logger.getLogger(AbstractCylindricalProjection.class.getName());
    
    /**
     * AbstractProjection name.
     */
    public static final String NAME = "Cylindrical projections";
    
    /**
     * Native longitude value in radians for cylindrical projection.
     */
    public static final double DEFAULT_PHI0 = 0;
    /**
     * Native latitude value in radians for cylindrical projection.
     */
    public static final double DEFAULT_THETA0 = 0;
    /**
     * Native longitude in radians of the ﬁducial point for the Cylindrical
     * projection.
     */
    private double phi0;
    /**
     * Native latitude in radians of the ﬁducial point for the Cylindrical 
     * projection.
     */
    private double theta0;

    /**
     * Creates a new cylindrical projection.
     * 
     * @param crval1 Celestial longitude in degrees of the ﬁducial point
     * @param crval2 Celestial latitude in degrees of the ﬁducial point
     */
    protected AbstractCylindricalProjection(final double crval1, final double crval2) {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[deg] (crval1,crval2) = ({0},{1})", new Object[]{crval1, crval2});
        setPhi0(DEFAULT_PHI0);
        setTheta0(DEFAULT_THETA0);
        setPhip(computeDefaultValueForPhip());        
        LOG.log(Level.FINEST, "(phi0,theta0)[DEG]=({0},{1})", new Object[]{FastMath.toDegrees(DEFAULT_PHI0), FastMath.toDegrees(DEFAULT_THETA0)});
        LOG.log(Level.FINEST, "phip[deg]={0}", FastMath.toDegrees(computeDefaultValueForPhip()));        
    }
    
    @Override
    public String getNameFamily() {
        return NAME;
    }  

    @Override
    public final void setPhi0(final double phi0) {
        this.phi0 = phi0;
    } 
    
    @Override
    public double getPhi0() {
        return phi0;
    }
    
    @Override
    public final void setTheta0(final double theta0) {
        this.theta0 = theta0;
    }     

    @Override
    public double getTheta0() {
        return theta0;
    }
    
    @Override
    public boolean inside(final double lon, final double lat) {      
       return true;
    }     
    
    @Override
    public boolean isLineToDraw(final double[] pos1, final double[] pos2) {
        return FastMath.abs(pos1[0] - pos2[0]) < 100;
    }    
    
    @Override
    public ProjectionParameter[] getProjectionParameters() {
        return new ProjectionParameter[]{};
    }
    
    @Override
    public final Logger getLogger() {
        return LOG;
    }    

}
