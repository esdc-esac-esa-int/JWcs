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
import io.github.malapert.jwcs.datum.CoordinateReferenceFrame;
import io.github.malapert.jwcs.proj.exception.JWcsError;
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.createRealIdentityMatrix;
import org.apache.commons.math3.linear.RealMatrix;


/**
 * Supergalactic coordinates are coordinates in a spherical coordinate system
 * which was designed to have its equator aligned with the supergalactic plane,
 * a major structure in the local universe formed by the preferential
 * distribution of nearby galaxy clusters (such as the Virgo cluster, 
 * the Great Attractor and the Pisces-Perseus supercluster) towards a 
 * (two-dimensional) plane. 
 * 
 * <p>The supergalactic plane was recognized by Gérard de Vaucouleurs in 1953 
 * from the Shapley-Ames Catalog, although a flattened distribution of nebulae
 * had been noted by William Herschel over 200 years earlier.
 * By convention, supergalactic latitude and supergalactic longitude are usually
 * denoted by SGB and SGL, respectively, by analogy to b and l conventionally
 * used for galactic coordinates. The zero point for supergalactic longitude
 * is defined by the intersection of this plane with the galactic plane.
 * 
 * <p>By convention, supergalactic latitude and supergalactic longitude are 
 * usually denoted by SGB and SGL, respectively, by analogy to b and l 
 * conventionally used for galactic coordinates. The zero point for 
 * supergalactic longitude is defined by the intersection of this plane 
 * with the galactic plane.
 *         
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 * @see <a href="https://en.wikipedia.org/wiki/Supergalactic_coordinate_system">Supergalactic coordinate system</a>
 */
public class SuperGalactic extends AbstractCrs {
    /**
     * Name of this coordinate system.
     */
    private static final CoordinateReferenceSystem SKY_NAME = CoordinateReferenceSystem.SUPER_GALACTIC;
    
    /**
     * The coordinate system.
     */
    private final AbstractCs coordinateSystem;
    
    public SuperGalactic() {
        this.coordinateSystem = CsFactory.create(AbstractCs.CoordinateSystem.SPHERICAL2D);
        this.coordinateSystem.getAxes()[0] = new AbstractCs.Axis("SGL", "Longitude super-galactic", AbstractCs.AxisDirection.EAST, AbstractCs.Unit.DEG);
        this.coordinateSystem.getAxes()[1] = new AbstractCs.Axis("SGB", "Latitude super-galactic", AbstractCs.AxisDirection.NORTH, AbstractCs.Unit.DEG);                       
    }
    
    /**
     * Returns the rotation matrix to convert from this current CRS to 
     * another one.
     * 
     * <p>The algorithm for the conversion is the following:
     * <ul>
     * <li>For {@link Equatorial} : <code>m = m3.multiply(m2).multiply(m1)</code> with <br>
     *  - <code>m1 = convertMatrixGal2Sgal().transpose()</code><br>
     *  - <code>m2 = convertMatrixEqB19502Gal().transpose()</code><br>
     *  - <code>m3 = convertMatrixEpoch12Epoch2(1950.0d, targetCrs.getEquinox(), CoordinateReferenceFrame.ReferenceFrame.FK4, targetCrs.getReferenceFrame(), Double.NaN)</code>
     * </li>
     * <li>For {@link Galactic} : <code>m = convertMatrixGal2Sgal().transpose()</code></li>
     * <li>For {@link SuperGalactic}: <code>m = createRealIdentityMatrix(3)</code></li>
     * <li>For {@link Ecliptic} : <code>m = m4.multiply(m3).multiply(m2).multiply(m1)</code> with <br>
     *  - <code>m1 = convertMatrixGal2Sgal().transpose()</code><br>
     *  - <code>m2 = convertMatrixEqB19502Gal().transpose()</code><br>
     *  - <code>m3 = convertMatrixEq2Ecl(targetCrs.getEquinox(), crs.getCoordinateReferenceFrame().getReferenceFrame())</code><br>
     *  - <code>m4 = convertMatrixEq2Ecl(targetCrs.getEquinox(), targetCrs.getReferenceFrame())</code>
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
     * @see NumericalUtility#createRealIdentityMatrix
     * @see AbstractCrs#convertMatrixGal2Sgal
     */    
    @Override
    protected RealMatrix getRotationMatrix(final AbstractCrs crs) throws JWcsError {
        final RealMatrix m;
        final CoordinateReferenceFrame targetCrs = crs.getCoordinateReferenceFrame();        
        final CoordinateReferenceSystem cs = crs.getCoordinateReferenceSystem();
        switch (cs) {
            case EQUATORIAL:
                RealMatrix m1 = convertMatrixGal2Sgal().transpose(); 
                RealMatrix m2 = convertMatrixEqB19502Gal().transpose();
                RealMatrix m3 = convertMatrixEpoch12Epoch2(1950.0d, targetCrs.getEquinox(), CoordinateReferenceFrame.ReferenceFrame.FK4, targetCrs.getReferenceFrame(), Double.NaN);
                m = m3.multiply(m2).multiply(m1);
                break;
            case GALACTIC:
                m = convertMatrixGal2Sgal().transpose();       
                break;
            case SUPER_GALACTIC:
                m = createRealIdentityMatrix(3);
                break;
            case ECLIPTIC:
                m1 = convertMatrixGal2Sgal().transpose();
                m2 = convertMatrixEqB19502Gal().transpose();
                m3 = convertMatrixEpoch12Epoch2(1950.0d, targetCrs.getEquinox(), CoordinateReferenceFrame.ReferenceFrame.FK4, targetCrs.getReferenceFrame(), Double.NaN);
                final RealMatrix m4 = convertMatrixEq2Ecl(targetCrs.getEquinox(), targetCrs.getReferenceFrame());
                m = m4.multiply(m3).multiply(m2).multiply(m1);
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
    public String toString() {
        return SKY_NAME.name();
    }

    @Override
    public CoordinateReferenceFrame getCoordinateReferenceFrame() {
        return null;
    }

    @Override
    public void setCoordinateReferenceFrame(final CoordinateReferenceFrame coordinateReferenceFrame) {
        // Intentionally empty because the coordinate is not time dependant.
    }    
}
