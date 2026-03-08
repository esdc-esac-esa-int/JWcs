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
package io.github.malapert.jwcs.crs;

import io.github.malapert.jwcs.coordsystem.AbstractCs;
import io.github.malapert.jwcs.coordsystem.AbstractCs.Axis;
import io.github.malapert.jwcs.coordsystem.CsFactory;
import io.github.malapert.jwcs.datum.FK4;
import io.github.malapert.jwcs.datum.CoordinateReferenceFrame;
import io.github.malapert.jwcs.datum.FK5;
import io.github.malapert.jwcs.proj.exception.JWcsError;
import java.util.Objects;
import org.apache.commons.math3.linear.RealMatrix;


/**
 * The ecliptic is the apparent path of the Sun on the celestial sphere, and is
 * the basis for the ecliptic coordinate system.
 *
 * <p>It also refers to the plane of this path, which is coplanar with both the
 * orbit of the Earth around the Sun and the apparent orbit of the Sun around
 * the Earth. The path of the Sun is not normally noticeable from the Earth's
 * surface because the Earth rotates, carrying the observer through the cycle of
 * sunrise and sunset, obscuring the apparent motion of the Sun with respect to
 * the stars.
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class Ecliptic extends AbstractCrs implements CoordinateReferenceFrame {

    /**
     * This coordinate system name.
     */
    private static final CoordinateReferenceSystem SKY_NAME = CoordinateReferenceSystem.ECLIPTIC;

    /**
     * The coordinate reference system.
     */
    private CoordinateReferenceFrame coordinateReferenceFrame;
    
    /**
     * The coordinate system.
     */
    private final AbstractCs coordinateSystem;

    /**
     * Creates the Ecliptic coordinate system based on the reference system.
     *
     * @param coordinateReferenceFrame the reference system
     */
    public Ecliptic(final CoordinateReferenceFrame coordinateReferenceFrame) {
        this.coordinateReferenceFrame = coordinateReferenceFrame;
        this.coordinateSystem = CsFactory.create(AbstractCs.CoordinateSystem.SPHERICAL2D);
        this.coordinateSystem.getAxes()[0] = new Axis("\u03BB", "Geocentric longitude", AbstractCs.AxisDirection.EAST, AbstractCs.Unit.DEG);
        this.coordinateSystem.getAxes()[1] = new Axis("\u03B2", "Geocentric latitude", AbstractCs.AxisDirection.NORTH, AbstractCs.Unit.DEG);
    }

    /**
     * Creates the Ecliptic coordinate system based on FK5 reference frame.
     */
    public Ecliptic() {
        this(new FK5());
    }

    /**
     * Transforms the coordinate reference frame if necessary. 
     * 
     * <p>In 'Representations of celestial coordinates in FITS' (Calabretta and
     * Greisen) we read that all reference systems are allowed for both
     * equatorial- and ecliptic coordinates, except FK4-NO-E, which is only
     * allowed for equatorial coordinates. If FK4-NO-E is given in combination
     * with an ecliptic sky system then silently FK4 is assumed.
     *
     * @param refSystem the coordinate reference system
     * @return the coordinate reference system
     */
    protected CoordinateReferenceFrame init(final CoordinateReferenceFrame refSystem) {
        CoordinateReferenceFrame result;
        if (ReferenceFrame.FK4_NO_E.equals(refSystem.getReferenceFrame())) {
            result = new FK4();
            result.setEpochObs(refSystem.getEpochObs());
            result.setEquinox(refSystem.getEquinox());
        } else {
            result = refSystem;
        }
        return result;
    }

    /**
     * Returns the rotation matrix to convert from this current CRS to 
     * another one.
     * 
     * <p>The algorithm for the conversion is the following:
     * <ul>
     * <li>For {@link Equatorial} : <code>m = m2.multiply(m1)</code> with <br>
     *  - <code>m1 = convertMatrixEq2Ecl(this.getEquinox(), getReferenceFrame()).transpose()</code><br>
     *  - <code>m2 = convertMatrixEpoch12Epoch2(this.getEquinox(), targetCrs.getEquinox(), getReferenceFrame(), targetCrs.getReferenceFrame(), Double.NaN);</code>
     * </li>
     * <li>For {@link Galactic} : <code>m = m3.multiply(m2).multiply(m1)</code> with <br>
     *  - <code>m1 = convertMatrixEq2Ecl(this.getEquinox(), getReferenceFrame()).transpose()</code><br>
     *  - <code>m2 = convertMatrixEpoch12Epoch2(this.getEquinox(), 1950.0d, getReferenceFrame(), CoordinateReferenceFrame.ReferenceFrame.FK4, Double.NaN)</code><br>
     *  - <code>m3 = convertMatrixEqB19502Gal()</code>
     * </li>
     * <li>For {@link SuperGalactic}: <code>m = m4.multiply(m3).multiply(m2).multiply(m1)</code> with <br>
     *  - <code>m1 = convertMatrixEq2Ecl(this.getEquinox(), getReferenceFrame()).transpose()</code><br>
     *  - <code>m2 = convertMatrixEpoch12Epoch2(this.getEquinox(), 1950.0d, getReferenceFrame(), CoordinateReferenceFrame.ReferenceFrame.FK4, Double.NaN)</code><br>
     *  - <code>m3 = convertMatrixEqB19502Gal()</code><br>
     *  - <code>m4 = convertMatrixGal2Sgal()</code>     
     * </li>
     * <li>For {@link Ecliptic} : <code>m = m3.multiply(m2).multiply(m1)</code> with <br>
     *  - <code>m1 = convertMatrixEq2Ecl(this.getEquinox(), getReferenceFrame()).transpose()</code><br>
     *  - <code>m2 = convertMatrixEpoch12Epoch2(this.getEquinox(), targetCrs.getEquinox(), getReferenceFrame(), targetCrs.getReferenceFrame(), Double.NaN)</code><br>
     *  - <code>m3 = convertMatrixEq2Ecl(this.getEquinox(), targetCrs.getReferenceFrame())</code>
     * </li>
     * </ul>     
     * 
     * @param crs target coordinate reference system
     * @return the rotation matrix to convert from this current CRS to 
     * another one
     * @throws JWcsError Unknown output crs
     * @see AbstractCrs#convertMatrixEq2Ecl
     * @see AbstractCrs#convertMatrixEpoch12Epoch2
     * @see AbstractCrs#convertMatrixEqB19502Gal      
     */
    @Override
    protected RealMatrix getRotationMatrix(final AbstractCrs crs) throws JWcsError {
        final RealMatrix m;
        final CoordinateReferenceFrame targetCrs = crs.getCoordinateReferenceFrame();
        final CoordinateReferenceSystem cs = crs.getCoordinateReferenceSystem();
        switch(cs) {
            case EQUATORIAL:
                RealMatrix m1 = convertMatrixEq2Ecl(this.getEquinox(), getReferenceFrame()).transpose();
                RealMatrix m2 = convertMatrixEpoch12Epoch2(this.getEquinox(), targetCrs.getEquinox(), getReferenceFrame(), targetCrs.getReferenceFrame(), Double.NaN);
                m = m2.multiply(m1);                
                break;
            case GALACTIC:
                m1 = convertMatrixEq2Ecl(this.getEquinox(), getReferenceFrame()).transpose();
                m2 = convertMatrixEpoch12Epoch2(this.getEquinox(), 1950.0d, getReferenceFrame(), CoordinateReferenceFrame.ReferenceFrame.FK4, Double.NaN);
                RealMatrix m3 = convertMatrixEqB19502Gal();
                m = m3.multiply(m2).multiply(m1);                
                break;
            case SUPER_GALACTIC:
                m1 = convertMatrixEq2Ecl(this.getEquinox(), getReferenceFrame()).transpose();
                m2 = convertMatrixEpoch12Epoch2(this.getEquinox(), 1950.0d, getReferenceFrame(), CoordinateReferenceFrame.ReferenceFrame.FK4, Double.NaN);
                m3 = convertMatrixEqB19502Gal();
                final RealMatrix m4 = convertMatrixGal2Sgal(); 
                m = m4.multiply(m3).multiply(m2).multiply(m1);
                break;
            case ECLIPTIC:
                m1 = convertMatrixEq2Ecl(this.getEquinox(), getReferenceFrame()).transpose();
                m2 = convertMatrixEpoch12Epoch2(this.getEquinox(), targetCrs.getEquinox(), getReferenceFrame(), targetCrs.getReferenceFrame(), Double.NaN);
                m3 = convertMatrixEq2Ecl(this.getEquinox(), targetCrs.getReferenceFrame());
                m = m3.multiply(m2).multiply(m1);                
                break;
            default:
                throw new JWcsError(String.format("Unknown output crs: %s", crs.getCoordinateReferenceSystem()));                                
        }
        return m;
    }

    @Override
    public CoordinateReferenceSystem getCoordinateReferenceSystem() {
        return SKY_NAME;
    }

    @Override
    public double getEquinox() {
        return this.getCoordinateReferenceFrame().getEquinox();
    }

    @Override
    public double getEpochObs() {
        return this.getCoordinateReferenceFrame().getEpochObs();
    }

    @Override
    public void setEquinox(final String equinox) {
        this.getCoordinateReferenceFrame().setEquinox(equinox);
    }

    @Override
    public void setEpochObs(final String epoch) {
        this.getCoordinateReferenceFrame().setEpochObs(epoch);
    }

    @Override
    public void setEquinox(final double equinox) {
        this.getCoordinateReferenceFrame().setEquinox(equinox);
    }

    @Override
    public void setEpochObs(final double epoch) {
        this.getCoordinateReferenceFrame().setEpochObs(epoch);
    }

    @Override
    public CoordinateReferenceFrame.ReferenceFrame getReferenceFrame() {
        return this.getCoordinateReferenceFrame().getReferenceFrame();
    }

    /**
     * Returns the coordinate reference frame.
     *
     * @return the coordinate reference frame
     */
    @Override
    public CoordinateReferenceFrame getCoordinateReferenceFrame() {
        return coordinateReferenceFrame;
    }

    /**
     * Sets the reference system.
     *
     * @param coordinateReferenceFrame the coordinateReferenceFrame to set
     */
    @Override
    public void setCoordinateReferenceFrame(final CoordinateReferenceFrame coordinateReferenceFrame) {
        this.coordinateReferenceFrame = coordinateReferenceFrame;
    }

    @Override
    public String toString() {
        return SKY_NAME.name();
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 11 * hash + Objects.hashCode(this.coordinateReferenceFrame);
        return hash;
    }

    @Override
    public boolean equals(final Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Ecliptic other = (Ecliptic) obj;
        return Objects.equals(this.coordinateReferenceFrame, other.coordinateReferenceFrame);          
    }

}
