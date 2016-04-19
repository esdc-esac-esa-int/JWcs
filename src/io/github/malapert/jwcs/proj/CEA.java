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

import io.github.malapert.jwcs.proj.exception.BadProjectionParameterException;
import io.github.malapert.jwcs.proj.exception.ProjectionException;
import io.github.malapert.jwcs.utility.NumericalUtils;

/**
 * The cylindrical equal area projection.
 * 
 * <p>
 * Reference: "Representations of celestial coordinates in FITS", 
 * M. R. Calabretta and E. W. Greisen
 * </p>
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class CEA extends CylindricalProjection {

    /**
     * Projection's name.
     */
    private static final String NAME_PROJECTION = "Cylindrical equal area";

    /**
     * Projection's description.
     */
    private static final String DESCRIPTION = "\u03BB=%s";

    /**
     * \u03BB Scaling parameter.
     */
    private double lambda;

    /**
     * Default value for \u03BB.
     */
    private static final int DEFAULT_VALUE = 1;

    /**
     * Constructs a CEA based on the celestial longitude and latitude of the
     * fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     *
     * \u03BB is set to {@link CEA#DEFAULT_VALUE}.
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @throws BadProjectionParameterException
     */
    public CEA(double crval1, double crval2) throws BadProjectionParameterException {
        this(crval1, crval2, DEFAULT_VALUE);
    }

    /**
     * Constructs a CEA based on the celestial longitude and latitude of the
     * fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) and \u03BB.
     * 
     * \u03BB is set by the FITS keyword PV<code>nbAxis</code>_1.
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @param lambda \u03BB dimensionless.
     * @throws BadProjectionParameterException
     */
    public CEA(double crval1, double crval2, double lambda) throws BadProjectionParameterException {
        super(crval1, crval2);
        if (NumericalUtils.equal(lambda, 0, DOUBLE_TOLERANCE) || lambda < 0 || lambda > 1.0) {
            throw new BadProjectionParameterException("CEA: Bad projection parameter for lambda " + lambda + " - lambda outside of range (0,1]");
        }
        setLambda(lambda);
    }

    @Override
    protected double[] project(double x, double y) throws ProjectionException {
        double xr = Math.toRadians(x);
        double yr = Math.toRadians(y);
        double phi = xr;
        double arg = getLambda() * yr;
        double theta = NumericalUtils.aasin(arg);
        double[] pos = {phi, theta};
        return pos;
    }

    @Override
    protected double[] projectInverse(double phi, double theta) throws ProjectionException {
        phi = phiRange(phi);
        double x = Math.toDegrees(phi);
        double y = Math.toDegrees(Math.sin(theta) / getLambda());
        double[] coord = {x, y};
        return coord;
    }

    /**
     * Returns \u03BB, the scaling parameter.
     *
     * @return the lambda
     */
    private double getLambda() {
        return lambda;
    }

    /**
     * Sets \u03BB, the scaling parameter.
     *
     * @param lambda the lambda to set
     */
    private void setLambda(double lambda) {
        this.lambda = lambda;
    }

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return String.format(DESCRIPTION, NumericalUtils.round(this.lambda));
    }

}
