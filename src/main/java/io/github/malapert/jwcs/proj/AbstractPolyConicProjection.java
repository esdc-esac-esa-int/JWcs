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

import io.github.malapert.jwcs.utility.NumericalUtility;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.util.FastMath;

/**
 * Polyconic projection.
 *
 * <p>Polyconics are generalizations of the standard conic projections; the
 * parallels of latitude are projected as circular arcs which may or may not be
 * concentric, and meridians are curved rather than straight as in the standard
 * conics.
 * 
 * <p>Ref : "Representations of celestial coordinates in FITS", Calabretta, M.R.,
 * and Greisen, E.W., (2002), Astronomy and Astrophysics, 395, 1077-1122. - p21
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public abstract class AbstractPolyConicProjection extends AbstractProjection {

    /**
     * Logger.
     */
    protected static final Logger LOG = Logger.getLogger(AbstractPolyConicProjection.class.getName());

    /**
     * AbstractProjection name.
     */
    public static final String NAME = "Polyconic and pseudoconic projections";

    /**
     * Native longitude value in radians for cylindrical projection.
     */
    protected static final double DEFAULT_PHI0 = 0;
    /**
     * Native latitude value in radians for cylindrical projection.
     */
    protected static final double DEFAULT_THETA0 = 0;

    /**
     * theta1.
     */
    private final double theta1;
    /**
     * Native longitude in radians of the ﬁducial point for the polyconic projection.
     */
    private double phio;
    /**
     * Native latitude in radians of the ﬁducial point for the polyconic projection.
     */
    private double theta0;

    /**
     * Creates a new polyconic projection.
     *
     * @param crval1 Celestial longitude in degrees of the ﬁducial point
     * @param crval2 Celestial latitude in degrees of the ﬁducial point
     * @param theta1 theta1 in degrees
     */
    protected AbstractPolyConicProjection(final double crval1, final double crval2, final double theta1) {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[deg] (crval1,crval2,theta1) = ({0},{1},{2})", new Object[]{crval1, crval2, theta1});

        if (NumericalUtility.equalValues(theta1, 0)) {
            this.theta1 = FastMath.toRadians(45);
            LOG.log(Level.WARNING,"theta1=0 not allowed, reseting to 45");
        } else {
            this.theta1 = FastMath.toRadians(theta1);
        }
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
    public double getPhi0() {
        return phio;
    }

    @Override
    public double getTheta0() {
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
     * Returns theta1 in radians.
     *
     * @return the theta1
     */
    protected double getTheta1() {
        return theta1;
    }

    @Override
    public boolean inside(final double lon, final double lat) {
        LOG.log(Level.FINER, "true");
        return true;
    }

    @Override
    public boolean isLineToDraw(final double[] pos1, final double[] pos2) {
        LOG.log(Level.FINER, "(pos1,pos2)=({0},{1}) ({2},{3})", new Object[]{pos1[0],pos1[1],pos2[0],pos2[1]});
        return FastMath.abs(pos1[0] - pos2[0]) < 50;
    }
    
    @Override
    public final Logger getLogger() {
        return LOG;
    }    
}
