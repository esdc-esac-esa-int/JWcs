/* 
 * Copyright (C) 2014 Jean-Christophe Malapert
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

import io.github.malapert.jwcs.proj.exception.BadProjectionParameterException;
import io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException;
import io.github.malapert.jwcs.utility.NumericalUtils;

/**
 * Conic equal area.
 *
 * <p>
 * The standard parallels in Alber's conic equal area projection are projected
 * as concentric arcs at their true length and separated so that the area
 * between them is the same as the corresponding area on the sphere. The other
 * parallels are then drawn as concentric arcs spaced so as to preserve the
 * area.
 * </p>
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 1.0
 */
public class COE extends ConicProjection {

    /**
     * Projection's name.
     */
    private static final String NAME_PROJECTION = "Conic equal area";
    
    /**
     * Projection's description.
     */
    private static final String DESCRIPTION = "\u03B8a=%s \u03B7=%s"; 
    
    /**
     * Creates an instance.
     *
     * @param crval1 Celestial longitude in degrees of the ﬁducial point
     * @param crval2 Celestial latitude in degrees of the ﬁducial point
     * @param theta_a (theta1 + theta2) / 2 in degrees
     * @param eta abs(theta1 - theta2) / 2 in degrees
     */
    public COE(double crval1, double crval2, double theta_a, double eta) {
        super(crval1, crval2, theta_a, eta);
    }

    /**
     * Computes the native spherical coordinates from the projection plane
     * coordinates.
     *
     *
     * @param x projection plane coordinate along X
     * @param y projection plane coordinate along Y
     * @return the native spherical coordinates in radians
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException when a projection parameter is wrong
     * @throws io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException When the pixel is beyond the visible projection
     */
    @Override
    protected double[] project(double x, double y) throws BadProjectionParameterException, PixelBeyondProjectionException {
        double xr = Math.toRadians(x);
        double yr = Math.toRadians(y);
        double theta1 = getTheta_a() - getEta();
        double theta2 = getTheta_a() + getEta();
        double gamma = Math.sin(theta1) + Math.sin(theta2);
        if (gamma == 0) {
            throw new BadProjectionParameterException("Projection parameters: sin(theta1) + sin(theta2) = 0");
        }
        double c = gamma * 0.5;
        
        
        double y0 = Math.sqrt(1.0d + Math.sin(theta1) * Math.sin(theta2) - gamma * Math.sin((theta1+theta2)*0.5)) / c;
        
        int sign = (getTheta_a() < 0) ? -1 : 1;
        double r_theta = sign * Math.sqrt(Math.pow(xr, 2) + Math.pow((y0 - yr), 2));
        double phi;
        if (NumericalUtils.equal(r_theta, 0, DOUBLE_TOLERANCE)) {
            phi = 0;
        } else {
            phi = NumericalUtils.aatan2(xr / r_theta, (y0 - yr) / r_theta) / c;
        }

        double w = 1.0d / gamma + Math.sin(theta1) * Math.sin(theta2) / gamma - gamma * Math.pow(r_theta * 0.5, 2);
        double theta;
        if (Math.abs(w) > 1) {
            if (NumericalUtils.equal(w, 1, DOUBLE_TOLERANCE)) {
                theta = HALF_PI;
            } else if (NumericalUtils.equal(w, -1, DOUBLE_TOLERANCE)) {
                theta = -HALF_PI;
            } else {
                throw new PixelBeyondProjectionException(
		 "COE: Calculation failed for (x,y) = (" + x + ", " + y + ")");
            }
        } else {
            theta = NumericalUtils.aasin(w);
        }
        double[] pos = {phi, theta};
        return pos;
    }

    /**
     * Computes the projection plane coordinates from the native spherical
     * coordinates.
     *
     * @param phi native spherical coordinate in radians along longitude
     * @param theta native spherical coordinate in radians along latitude
     * @return the projection plane coordinates
     */
    @Override
    protected double[] projectInverse(double phi, double theta) {
        double theta1 = getTheta_a() - getEta();
        double theta2 = getTheta_a() + getEta();
        double gamma = Math.sin(theta1) + Math.sin(theta2);
        double c = gamma * 0.5;
        double y0 = Math.sqrt(1 + Math.sin(theta1) * Math.sin(theta2) - gamma * Math.sin(getTheta_a())) / c;
        phi = phi % (Math.PI * 2);
        if (phi > Math.PI) {
            phi -= Math.PI * 2;
        } else if (phi < -Math.PI) {
            phi += Math.PI * 2;
        }
        double r_theta = Math.sqrt(1.0d + Math.sin(theta1) * Math.sin(theta2) - gamma * Math.sin(theta)) / c;
        double x = Math.toDegrees(r_theta * Math.sin(c * phi));
        double y = Math.toDegrees(-r_theta * Math.cos(c * phi) + y0);

        double[] coord = {x, y};
        return coord;
    }
    
    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return String.format(DESCRIPTION, NumericalUtils.round(Math.toDegrees(this.getTheta_a())), NumericalUtils.round(Math.toDegrees(this.getEta())));
    }    

}
