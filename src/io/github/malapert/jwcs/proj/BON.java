/* 
 * Copyright (C) 2014-2016 Jean-Christophe Malapert
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
public class BON extends AbstractPolyConicProjection {
    
    /**
     * Projection's name.
     */
    private final static String NAME_PROJECTION = "Bonne’s equal area";
    
    /**
     * Projection's description.
     */
    private final static String DESCRIPTION = "no limits";     

    /**
     * SFL projection.
     */
    private SFL sfl;

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
        if (NumericalUtility.equal(theta1,0)) {
            this.sfl = new SFL(crval1, crval2);
        }
    }

    @Override
    protected double[] project(final double x, final double y) {
        LOG.log(Level.FINER, "INPUTS[Deg] (x,y)=({0},{1})", new Object[]{x,y});                        
        final double[] result;
        if (this.sfl == null) {
            final double xr = FastMath.toRadians(x);
            final double yr = FastMath.toRadians(y);
            final double y0 = getTheta1() + 1.0d / FastMath.tan(getTheta1());            
            final double r_theta = FastMath.signum(getTheta1())* FastMath.sqrt(FastMath.pow(xr, 2) + FastMath.pow(y0 - yr, 2));
            final double aphi = NumericalUtility.aatan2(xr / r_theta, (y0 - yr) / r_theta);
            
            final double theta = y0 - r_theta;
            final double cos_theta = FastMath.cos(theta);
            final double phi;
            if (NumericalUtility.equal(cos_theta,0)) {
                phi = 0;
            } else {
                phi = aphi * r_theta / cos_theta;
            }
            final double[] pos = {phi, theta};
            result = pos;
        } else {
            result = this.sfl.project(x, y);
        }
        LOG.log(Level.FINER, "OUTPUTS[Deg] (phi,theta)=({0},{1})", new Object[]{FastMath.toDegrees(result[0]),FastMath.toDegrees(result[1])});                                
        return result;
    }

    @Override
    protected double[] projectInverse(final double phi, final double theta) {
        LOG.log(Level.FINER, "INPUTS[Deg] (phi,theta)=({0},{1})", new Object[]{FastMath.toDegrees(phi),FastMath.toDegrees(theta)});                                        
        double[] result;
        if (sfl == null) {
            final double y0 = getTheta1() + 1.0d / FastMath.tan(getTheta1());
            final double r_theta = y0 - theta;
            final double aphi;
            if (NumericalUtility.equal(r_theta, 0)) {
                aphi=0;
            } else {
                aphi = phi * FastMath.cos(theta) / r_theta;
            }
            final double x = r_theta * FastMath.sin(aphi);
            final double y = -r_theta * FastMath.cos(aphi) + y0;
            final double[] coord = {FastMath.toDegrees(x), FastMath.toDegrees(y)};
            result = coord;
        } else {
            result = sfl.projectInverse(phi, theta);
        }
        LOG.log(Level.FINER, "OUTPUTS[Deg] (x,y)=({0},{1})", new Object[]{result[0],result[1]});                                
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
