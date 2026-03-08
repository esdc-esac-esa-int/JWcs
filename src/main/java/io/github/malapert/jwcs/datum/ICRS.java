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
package io.github.malapert.jwcs.datum;

import static io.github.malapert.jwcs.utility.TimeUtility.epochs;

/**
 * The International Celestial Reference System, for optical data 
 * realized through the Hipparcos catalog. 
 * 
 * <p>By definition, ICRS is not an equatorial system, but it is 
 * very close to the FK5 (J2000) system. No Equinox value is required.
 * 
 * <p>The International Celestial Reference System (ICRS) is the current standard
 * celestial reference system adopted by the International Astronomical Union
 * (IAU). Its origin is at the barycenter of the solar system, with axes that
 * are intended to be "fixed" with respect to space. ICRS coordinates are
 * approximately the same as equatorial coordinates: the mean pole at J2000.0
 * in the ICRS lies at 17.3±0.2 mas in the direction 12 h and 5.1±0.2 mas
 * in the direction 18 h. The mean equinox of J2000.0 is shifted from
 * the ICRS right ascension origin by 78±10 mas (direct rotation around
 * the polar axis).
 *         
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class ICRS implements CoordinateReferenceFrame {
  
    /**
     * The name of this reference frame.
     */
    private static final CoordinateReferenceFrame.ReferenceFrame REF_SYSTEM = CoordinateReferenceFrame.ReferenceFrame.ICRS;

    /**
     * The default value of the epoch sets to J2000.
     */
    private static final String DEFAULT_EPOCH = "J2000";
    
    /**
     * The epoch of the equinox.
     * 
     * <p>An equinox is an astronomical event in which the plane of Earth's 
     * equator passes through the center of the Sun, which occurs twice each 
     * year, around 20 March and 23 September.
     */
    private double equinox;

    /**
     * Creates a new ICRS reference frame based on {@link ICRS#DEFAULT_EPOCH}.
     */
    public ICRS() {
        init(DEFAULT_EPOCH);
    }

    /**
     * Initialization.
     * @param equinox the epoch of equinox
     */
    private void init(final String equinox) {
        this.equinox = epochs(equinox)[1];
    }

    @Override
    public CoordinateReferenceFrame.ReferenceFrame getReferenceFrame() {
        return REF_SYSTEM;
    }

    /**
     * Returns NaN.
     * 
     * <p>No need to specify an epoch of observation in ICRS reference frame
     */    
    @Override
    public double getEpochObs() {
        return Double.NaN;
    }

    /**
     * Returns the equinox as a Julian epoch.
     */    
    @Override
    public double getEquinox() {
        return this.equinox;
    }
    
    /**
     * Do nothing.
     * No need to specify the epoch of equinox
     */ 
    @Override   
    public void setEquinox(final String equinox) {
        //Intentionally empty because the equinox is invariant in ICRS frame.
    }    
    
    /**
     * Do nothing.
     * No need to specify the epoch of equinox
     */    
    @Override
    public void setEquinox(final double equinox) {
        //Intentionally empty because the equinox is invariant in ICRS frame.
    }      
    
    /**
     * Do nothing.
     * No need to specify the epoch of observation in FK5 reference frame
     */    
    @Override
    public void setEpochObs(final String epochObs) { 
        //Intentionally empty because the equinox is invariant in ICRS frame.
    }

    /**
     * Do nothing.
     * No need to specify the epoch of observation in FK5 reference frame
     */    
    @Override    
    public void setEpochObs(final double epochObs) {
        //Intentionally empty because the equinox is invariant in ICRS frame.
    }     

    @Override
    public String toString() {
        return "ICRS";
    } 

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 53 * hash + (int) (Double.doubleToLongBits(this.equinox) ^ (Double.doubleToLongBits(this.equinox) >>> 32));
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
        final ICRS other = (ICRS) obj;
        return  Double.doubleToLongBits(this.equinox) == Double.doubleToLongBits(other.equinox);
    }
    
}
