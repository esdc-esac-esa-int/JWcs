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

import io.github.malapert.jwcs.utility.NumericalUtils;

/**
 * Zenithal equidistant.
 *
 * <p>
 * The zenithal equidistant projection ﬁrst appeared in Greisen (1983) as ARC.
 * It is widely used as the approximate projection of Schmidt telescopes.
 * The native meridians are uniformly divided to give equispaced parallels.
 *  This projection was also known in antiquity.
 * 
 * @see <a href="http://www.atnf.csiro.au/people/mcalabre/WCS/ccs.pdf">
 * "Representations of celestial coordinates in FITS"</a>, 
 * M. R. Calabretta and E. W. Greisen - page 14
 * </p>
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class ARC extends ZenithalProjection {
    
    /**
     * Projection's name.
     */
    private static final String NAME_PROJECTION = "Zenithal equidistant";
    
    /**
     * Projection's description.
     */
    private static final String DESCRIPTION = "no limits";

    /**
     * Creates a new ARC projection based on the celestial longitude and 
     * latitude of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>)
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the fiducial point
     */    
    public ARC(double crval1, double crval2) {
        super(crval1, crval2);
    }
 
    @Override
    protected double[] project(double x, double y) {
        double xr = Math.toRadians(x);
        double yr = Math.toRadians(y);
        double r_theta = Math.sqrt(Math.pow(xr, 2) + Math.pow(yr, 2));
        double phi;
        if (NumericalUtils.equal(r_theta, 0, DOUBLE_TOLERANCE)) {
            phi = 0.0;
        } else {
            phi = NumericalUtils.aatan2(xr, -yr);
        }
        double theta = HALF_PI - r_theta;
        double[] pos = {phi, theta};
        return pos;
    }

    @Override
    protected double[] projectInverse(double phi, double theta) {
        phi = phiRange(phi);
        double r = HALF_PI - theta;  // theta between [-HALF_PI, HALF_PI] => no need to test
        double x = r * Math.sin(phi);
        double y = -r * Math.cos(phi);
        double[] pos = {Math.toDegrees(x), Math.toDegrees(y)};
        return pos;
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
