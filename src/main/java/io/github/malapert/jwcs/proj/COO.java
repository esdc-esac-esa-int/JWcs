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
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.HALF_PI;
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * Conic orthomorphic.
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class COO extends AbstractConicProjection {

    /**
     * Projection's name.
     */
    private static final String NAME_PROJECTION = " Conic orthomorphic";

    /**
     * Projection's description.
     */
    private static final String DESCRIPTION = "\u03B8a=%s \u03B7=%s";

    /**
     * Constant of the cone.
     */
    private final double c;

    /**
     * tan((HALF_PI - \u03B8<sub>1</sub>) * 0.5).
     */
    private final double tan1;

    /**
     * tan((HALF_PI - \u03B8<sub>2</sub>) * 0.5).
     */
    private final double tan2;
    
    /**
     * Constructs a COO projection based on the default celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) and 
     * 03B8<sub>a</sub> = 45 and \u03B7 = 25.
     *
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException When projection parameters are wrong
     */     
    public COO() throws BadProjectionParameterException {
        this(FastMath.toDegrees(AbstractConicProjection.DEFAULT_PHI0), 45, 45, 25);        
    }

    /**
     * Constructs a COO projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>) and
     * 03B8<sub>a</sub> and \u03B7.
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @param theta_a \u03B8<sub>a</sub> in degrees and defined as
     * \u03B8<sub>a</sub>=(\u03B8<sub>1</sub>+\u03B8<sub>2</sub>)/2
     * @param eta \u03B7 in degrees and defined as
     * \u03B7=|\u03B8<sub>1</sub>-\u03B8<sub>2</sub>|/2
     * @throws
     * io.github.malapert.jwcs.proj.exception.BadProjectionParameterException
     * When projection parameters are wrong
     */
    public COO(final double crval1, final double crval2, final double theta_a, final double eta) throws BadProjectionParameterException {
        super(crval1, crval2, theta_a, eta);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2,theta_a,eta)=({0},{1},{2},{3})", new Object[]{crval1, crval2, theta_a, eta});
        this.tan1 = FastMath.tan((HALF_PI - this.getTheta1()) * 0.5);
        this.tan2 = FastMath.tan((HALF_PI - this.getTheta2()) * 0.5);
        this.c = computeC(this.tan1, this.tan2);
    }

    /**
     * Computes the constant of the cone.
     *
     * <p>c = log(cos\u03B8<sub>2</sub> / cos\u03B8<sub>1</sub>) / log
     * ({@link COO#tan1} / {@link COO#tan2})
     *
     * @param tan1 {@link COO#tan1}
     * @param tan2 {@link COO#tan2}
     * @return the value
     * @throws BadProjectionParameterException (theta1,theta2). c must be != 0
     */
    private double computeC(final double tan1, final double tan2) throws BadProjectionParameterException {
        final double cte = NumericalUtility.equal(getTheta1(), getTheta2()) ? FastMath.sin(getTheta1()) : FastMath.log(FastMath.cos(getTheta2()) / FastMath.cos(getTheta1())) / FastMath.log(tan2 / tan1);
        if (NumericalUtility.equal(cte, 0)) {
            throw new BadProjectionParameterException(this, "(theta1,theta2). c must be != 0");
        }
        return cte;
    }

    @Override
    protected double[] project(final double x, final double y) {
        final double xr = FastMath.toRadians(x);
        final double yr = FastMath.toRadians(y);
        final double psi = NumericalUtility.equal(getTan1(), 0) ? FastMath.cos(getTheta2()) / (getC() * FastMath.pow(getTan2(), getC())) : FastMath.cos(getTheta1()) / (getC() * FastMath.pow(getTan1(), getC()));
        final double y0 = psi * FastMath.pow(FastMath.tan((HALF_PI - getThetaA()) * 0.5), getC());
        final double r_theta = FastMath.signum(getThetaA()) * FastMath.sqrt(FastMath.pow(xr, 2) + FastMath.pow(y0 - yr, 2));
        final double phi = computePhi(xr, yr, r_theta, y0, getC());
        final double theta = HALF_PI - 2 * FastMath.atan(FastMath.pow(r_theta / psi, 1.0 / getC()));
        return new double[] { phi, theta};
    }

    @Override
    protected double[] projectInverse(final double phi, final double theta) throws BadProjectionParameterException {
        final double psi = NumericalUtility.equal(getTan1(), 0) ? FastMath.cos(getTheta2()) / (getC() * FastMath.pow(getTan2(), getC())) : FastMath.cos(getTheta1()) / (getC() * FastMath.pow(getTan1(), getC()));
        if (NumericalUtility.equal(psi, 0)) {
            throw new BadProjectionParameterException(this, "(theta_a, eta) = (" + getThetaA() + ", " + getEta() + ")");
        }
        final double y0 = psi * FastMath.pow(FastMath.tan((HALF_PI - getThetaA()) * 0.5), getC());
        final double r_theta = psi * FastMath.pow(FastMath.tan((HALF_PI - theta) * 0.5), getC());
        final double x = computeX(phi, r_theta, getC());
        final double y = computeY(phi, r_theta, getC(), y0);
        return new double[] { FastMath.toDegrees(x), FastMath.toDegrees(y)};
    }

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return String.format(DESCRIPTION, NumericalUtility.round(FastMath.toDegrees(this.getThetaA())), NumericalUtility.round(FastMath.toDegrees(this.getEta())));
    }

    @Override
    public boolean inside(final double lon, final double lat) {
        LOG.log(Level.FINER, "(lon,lat)=({0},{1}) {2}", new Object[]{FastMath.toDegrees(lon), FastMath.toDegrees(lat), super.inside(lon, lat) && !NumericalUtility.equal(lat, -HALF_PI)});
        return super.inside(lon, lat) && !NumericalUtility.equal(lat, -HALF_PI);
    }

    @Override
    public ProjectionParameter[] getProjectionParameters() {
        final ProjectionParameter p1 = new ProjectionParameter("\u03B8a", AbstractJWcs.PV21, new double[]{-90, 90}, -45);
        final ProjectionParameter p2 = new ProjectionParameter("\u03B7", AbstractJWcs.PV22, new double[]{0, 90}, 0);
        return new ProjectionParameter[]{p1, p2};
    }

    /**
     * Returns c.
     *
     * @return the c
     */
    private double getC() {
        return c;
    }

    /**
     * Returns tan1.
     *
     * @return the tan1
     */
    private double getTan1() {
        return tan1;
    }

    /**
     * Returns tan2.
     *
     * @return the tan2
     */
    private double getTan2() {
        return tan2;
    }
}
