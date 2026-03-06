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
import org.apache.commons.math3.util.FastMath;

/**
 * The plate carrée projection.
 * 
 * <p>The equator and all meridians are correctly scaled in the plate
 * carrée projection.
 * 
 * <p>Reference: "Representations of celestial coordinates in FITS", 
 * M. R. Calabretta and E. W. Greisen - page 16
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public final class CAR extends AbstractCylindricalProjection{
    
    /**
     * Projection's name.
     */
    private final static String NAME_PROJECTION = "Plate carrée";
    
    /**
     * Projection's description.
     */
    private final static String DESCRIPTION = "no limits";
    
    /**
     * Constructs a CAR based on the default celestial longitude and 
     * latitude of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).     
     */    
    public CAR() {
        this(FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_THETA0));
    }

    /**
     * Constructs a CAR based on the celestial longitude and 
     * latitude of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     * 
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the fiducial point
     */
    public CAR(final double crval1, final double crval2) {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2)=({0},{1})", new Object[]{crval1,crval2});                
    }
    
    /**
     * Computes the native spherical coordinates (\u03D5, \u03B8) from the projection plane
     * coordinates (x, y).
     * 
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes \u03D5 = x</li>
     * <li>computes \u03B8 = y</li>
     * </ul>
     * 
     * @param x projection plane coordinate along X
     * @param y projection plane coordinate along Y
     * @return the native spherical coordinates (\u03D5, \u03B8) in radians
     */      
    @Override
    protected double[] project(final double x, final double y) {
        final double phi = FastMath.toRadians(x);
        final double theta = FastMath.toRadians(y);
        return new double[] { phi, theta};
    }

    /**
     * Computes the projection plane coordinates (x, y) from the native spherical
     * coordinates (\u03D5, \u03B8).
     *
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes x = \u03D5</li>
     * <li>computes y = \u03B8</li>
     * </ul>
     * 
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param theta the native spherical coordinate (\u03B8) in radians along latitude
     * @return the projection plane coordinates
     */    
    @Override
    protected double[] projectInverse(final double phi, final double theta) {
        final double x = FastMath.toDegrees(phi);
        final double y = FastMath.toDegrees(theta);
        return new double[] { x, y};
    }

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return DESCRIPTION;
    }

}
