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
import io.github.malapert.jwcs.utility.NumericalUtility;
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * Bonnes's equal area.
 *
 * <p>In Bonne's pseudoconic projection19 all parallels are projected as concentric
 * equidistant arcs of circles of true length and true spacing. This is sucient
 * to guarantee that it is an equal area projection.
 *
 * <p>Reference: "Representations of celestial coordinates in FITS", M. R.
 * Calabretta and E. W. Greisen - page 21
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public final class BON extends AbstractPolyConicProjection {
    
    /**
     * Projection's name.
     */
    private static final String NAME_PROJECTION = "Bonne’s equal area";
    
    /**
     * Projection's description.
     */
    private static final String DESCRIPTION = "no limits";

    /**
     * SFL projection.
     */
    private SFL sfl;

    /**
     * Constructs a BON projection by providing the default celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) and \u03B8<sub>1</sub>=45&deg;.          
     */    
    public BON() {
        this(AbstractPolyConicProjection.DEFAULT_PHI0, AbstractPolyConicProjection.DEFAULT_THETA0, 45.0);
    }

    /**
     * Constructs a BON projection by providing celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) and \u03B8<sub>1</sub>.
     *
     * <p>\u03B8<sub>1</sub> can be set by the FITS keyword PV<code>nbAxis</code>_1 in degrees.
     * 
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the fiducial point
     * @param theta1 \u03B8<sub>1</sub> latitude in degrees. 
     */
    public BON(final double crval1, final double crval2, final double theta1) {
        super(crval1, crval2, theta1);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2,theta1)=({0},{1},{2})", new Object[]{crval1,crval2,theta1});                
        if (NumericalUtility.equalValues(theta1, 0)) {
            this.sfl = new SFL(crval1, crval2);
        }
    }

    /**
     * Computes the native spherical coordinates (\u03D5, \u03B8) from the projection plane
     * coordinates (x, y).
     * 
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes radius : sign(\u03B8<sub>1</sub>) * sqrt(x<sup>2</sup> + (y0 - y)<sup>2</sup>) with y0 = {@link BON#computeY0() }</li>
     * <li>computes \u03B8 : y0 - radius</li>    
     * <li>computes \u03D5 : aphi * radius / cos\u03B8 with aphi = arg((y0 - y) / radius, x / radius)</li>
     * </ul>
     * 
     * <p>Special case : \u03B8<sub>1</sub>=0 is a special case, which is handled by {@link SFL}
     *
     * @param x projection plane coordinate along X
     * @param y projection plane coordinate along Y
     * @return the native spherical coordinates (\u03D5, \u03B8) in radians
     */     
    @Override
    protected double[] project(final double x, final double y) {
        final double[] result;
        if (this.sfl == null) {
            final double xr = FastMath.toRadians(x);
            final double yr = FastMath.toRadians(y);
            // compute radius
            final double y0 = computeY0();            
            final double r_theta = FastMath.signum(getTheta1())* FastMath.sqrt(FastMath.pow(xr, 2) + FastMath.pow(y0 - yr, 2));
            // compute theta
            final double theta = y0 - r_theta;
            // compute phi
            final double aphi = NumericalUtility.aatan2(xr / r_theta, (y0 - yr) / r_theta);
            final double cos_theta = FastMath.cos(theta);
            final double phi;
            if (NumericalUtility.equalValues(cos_theta, 0)) {
                phi = 0;
            } else {
                phi = aphi * r_theta / cos_theta;
            }
            result = new double[] { phi, theta};
        } else {
            result = this.sfl.project(x, y);
        }
        return result;
    }
    
    /**
     * Computes y0.
     * 
     * <p>y0 = \u03B8<sub>1</sub> + 1 / tan\u03B8<sub>1</sub>
     * 
     * @return y0
     */
    private double computeY0() {
        return getTheta1() + 1.0d / FastMath.tan(getTheta1());
    }

    /**
     * Computes the projection plane coordinates (x, y) from the native spherical
     * coordinates (\u03D5, \u03B8).
     *
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes radius : y0 - theta with y0 = {@link BON#computeY0() }</li>
     * <li>computes x : radius * sin(aphi)</li>
     * <li>computes y : -radius * cos(aphi) + y0 with aphi = phi * cos(\u03B8) / radius</li>
     * </ul>
     * 
     * <p>Special case : \u03B8<sub>1</sub>=0 is a special case, which is handled by {@link SFL}     
     * 
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param theta the native spherical coordinate (\u03B8) in radians along latitude
     * @return the projection plane coordinates
     */     
    @Override
    protected double[] projectInverse(final double phi, final double theta) {
        double[] result;
        if (sfl == null) {
            // computes radius
            final double y0 = computeY0();
            final double r_theta = y0 - theta;
            // computes x and y
            final double aphi;
            if (NumericalUtility.equalValues(r_theta, 0)) {
                aphi=0;
            } else {
                aphi = phi * FastMath.cos(theta) / r_theta;
            }
            final double x = r_theta * FastMath.sin(aphi);
            final double y = -r_theta * FastMath.cos(aphi) + y0;
            result = new double[] { FastMath.toDegrees(x), FastMath.toDegrees(y)};
        } else {
            result = sfl.projectInverse(phi, theta);
        }
        return result;
    }
    
    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return DESCRIPTION;
    }     

    @Override
    public ProjectionParameter[] getProjectionParameters() {
        final ProjectionParameter p1 = new ProjectionParameter("\u03B81", AbstractJWcs.PV21, new double[]{-90, 90}, 45);
        return new ProjectionParameter[]{p1};
    }

}
