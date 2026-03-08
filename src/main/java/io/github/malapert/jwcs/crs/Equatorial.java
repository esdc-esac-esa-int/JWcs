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
import io.github.malapert.jwcs.coordsystem.CsFactory;
import io.github.malapert.jwcs.datum.ICRS;
import io.github.malapert.jwcs.datum.CoordinateReferenceFrame;
import io.github.malapert.jwcs.proj.exception.JWcsError;
import java.util.Objects;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * The equatorial coordinate system is a widely used celestial coordinate system
 * used to specify the positions of celestial objects. 
 * 
 * <p>It may be implemented in spherical or rectangular coordinates, both defined
 * by an origin at the center of the Earth, a fundamental plane consisting of
 * the projection of the Earth's equator onto the celestial sphere
 * (forming the celestial equator), a primary direction towards the vernal
 * equinox, and a right-handed convention.
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class Equatorial extends AbstractCrs implements CoordinateReferenceFrame {

    /**
     * This coordinate system name.
     */
    private static final CoordinateReferenceSystem SKY_NAME = CoordinateReferenceSystem.EQUATORIAL;
    
    /**
     * The reference frame of the Equatorial coordinate system.
     */
    private CoordinateReferenceFrame coordinateReferenceFrame;

    /**
     * The coordinate system.
     */
    private final AbstractCs coordinateSystem;
    
    /**
     * Creates an Equatorial coordinate system based on the coordinate 
     * reference frame.
     * @param coordinateReferenceFrame the coordinate reference frame
     */    
    public Equatorial(final CoordinateReferenceFrame coordinateReferenceFrame) {
        this.coordinateReferenceFrame = coordinateReferenceFrame;
        this.coordinateSystem = CsFactory.create(AbstractCs.CoordinateSystem.SPHERICAL2D);
        this.coordinateSystem.getAxes()[0] = new AbstractCs.Axis("\u03B1", "Right ascension", AbstractCs.AxisDirection.EAST, AbstractCs.Unit.DEG);
        this.coordinateSystem.getAxes()[1] = new AbstractCs.Axis("\u1D5F", "Declination", AbstractCs.AxisDirection.NORTH, AbstractCs.Unit.DEG);        
    }

    /**
     * Creates an Equatorial coordinate system based on the ICRS reference frame.
     */
    public Equatorial() {
        this(new ICRS());
    }
   
    /**
     * Returns the rotation matrix to convert from this current CRS to 
     * another one.
     * 
     * <p>The algorithm for the conversion is the following:
     * <ul>
     * <li>For {@link Equatorial} : <code>m = convertMatrixEpoch12Epoch2(getEquinox(), targetCrs.getEquinox(), getReferenceFrame(), targetCrs.getReferenceFrame(), getEpochObs())</code></li>
     * <li>For {@link Galactic} : <code>m = m2.multiply(m1)</code> with : <br>
     *  - <code>m1 = convertMatrixEpoch12Epoch2(getEquinox(), 1950.0d, getReferenceFrame(), CoordinateReferenceFrame.ReferenceFrame.FK4, Double.NaN)</code><br>
     *  - <code>m2 = convertMatrixEqB19502Gal()</code>
     * </li>
     * <li>For {@link SuperGalactic}: <code>m3.multiply(m2).multiply(m1)</code> with <br>
     *  - <code>m1 = convertMatrixEpoch12Epoch2(getEquinox(), 1950.0d, getReferenceFrame(), CoordinateReferenceFrame.ReferenceFrame.FK4, Double.NaN)</code><br>
     *  - <code>m2 = convertMatrixEqB19502Gal()</code><br>
     *  - <code>m3 = convertMatrixGal2Sgal()</code>
     * </li>
     * <li>For {@link Ecliptic} : <code>m = m2.multiply(m1)</code> with <br>
     *  - <code>m1 = convertMatrixEpoch12Epoch2(getEquinox(), targetCrs.getEquinox(), getReferenceFrame(), targetCrs.getReferenceFrame(), Double.NaN)</code><br>
     *  - <code>m2 = convertMatrixEq2Ecl(targetCrs.getEquinox(), targetCrs.getReferenceFrame())</code>
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
                m = convertMatrixEpoch12Epoch2(getEquinox(), targetCrs.getEquinox(), getReferenceFrame(), targetCrs.getReferenceFrame(), getEpochObs());
                break;
            case GALACTIC:
                RealMatrix m1 = convertMatrixEpoch12Epoch2(getEquinox(), 1950.0d, getReferenceFrame(), CoordinateReferenceFrame.ReferenceFrame.FK4, Double.NaN);
                RealMatrix m2 = convertMatrixEqB19502Gal();
                m = m2.multiply(m1);                
                break;
            case SUPER_GALACTIC:
                m1 = convertMatrixEpoch12Epoch2(getEquinox(), 1950.0d, getReferenceFrame(), CoordinateReferenceFrame.ReferenceFrame.FK4, Double.NaN);
                m2 = convertMatrixEqB19502Gal();
                final RealMatrix m3 = convertMatrixGal2Sgal();
                m = m3.multiply(m2).multiply(m1);                
                break;
            case ECLIPTIC:
                m1 = convertMatrixEpoch12Epoch2(getEquinox(), targetCrs.getEquinox(), getReferenceFrame(), targetCrs.getReferenceFrame(), Double.NaN);
                m2 = convertMatrixEq2Ecl(targetCrs.getEquinox(), targetCrs.getReferenceFrame());
                m = m2.multiply(m1);                
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
     * @return the coordinate reference frame
     */
    @Override
    public CoordinateReferenceFrame getCoordinateReferenceFrame() {
        return coordinateReferenceFrame;
    }

    /**
     * Sets the coordinate reference system.
     * @param coordinateReferenceFrame the coordinateReferenceFrame to set
     */
    @Override
    public void setCoordinateReferenceFrame(final CoordinateReferenceFrame coordinateReferenceFrame) {
        this.coordinateReferenceFrame = coordinateReferenceFrame;
    }

    @Override
    public String toString() {
        return SKY_NAME+"("+coordinateReferenceFrame+")";
    } 

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 83 * hash + Objects.hashCode(this.coordinateReferenceFrame);
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
        final Equatorial other = (Equatorial) obj;
        return Objects.equals(this.coordinateReferenceFrame, other.coordinateReferenceFrame);
    }
    
}
