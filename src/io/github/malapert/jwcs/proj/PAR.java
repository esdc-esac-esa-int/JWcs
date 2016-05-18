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

import io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException;
import io.github.malapert.jwcs.utility.NumericalUtils;
import java.util.logging.Level;

/**
 * Parabolic.
 *
 * <p>The meridians are projected as parabolic arcs which intersect the poles and
 * correctly divide the equator, and the parallels of latitude are spaced so as
 * to make it an equal area projection. 
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class PAR extends CylindricalProjection {

    /**
     * Projection's name.
     */
    private final static String NAME_PROJECTION = "Parabolic";

    /**
     * Projection's description.
     */
    private final static String DESCRIPTION = "no limits";

    /**
     * Constructs a PAR projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     */
    public PAR(final double crval1, final double crval2) {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2)=({0},{1})", new Object[]{crval1,crval2});                                
    }

    @Override
    protected double[] project(final double x, final double y) throws PixelBeyondProjectionException {
        LOG.log(Level.FINER, "INPUTS[Deg] (x,y)=({0},{1})", new Object[]{x,y});                                                                                                
        final double xr = Math.toRadians(x);
        final double yr = Math.toRadians(y);
        final double theta = 3 * NumericalUtils.aasin(yr / Math.PI);
        if (Double.isNaN(theta)) {
            throw new PixelBeyondProjectionException(this, "(x,y)=(" + x + "," + y + ")");
        }
        final double phi = xr / (1 - 4 * Math.pow(yr / Math.PI, 2));
        final double[] pos = {phi, theta};
        LOG.log(Level.FINER, "OUTPUTS[Deg] (phi,theta)=({0},{1})", new Object[]{Math.toDegrees(phi),Math.toDegrees(theta)});                                                                                                        
        return pos;
    }

    @Override
    protected double[] projectInverse(final double phi, final double theta) {
        LOG.log(Level.FINER, "INPUTS[Deg] (phi,theta)=({0},{1})", new Object[]{Math.toDegrees(phi),Math.toDegrees(theta)});                                                                                                                
        final double y = Math.toDegrees(Math.PI * Math.sin(theta / 3d));
        final double x = Math.toDegrees(phi * (2d * Math.cos(theta / 1.5d) - 1d));
        final double[] coord = {x, y};
        LOG.log(Level.FINER, "OUTPUTS[Deg] (x,y)=({0},{1})", new Object[]{x,y});                                                                                                        
        return coord;
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
