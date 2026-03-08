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

import io.github.malapert.jwcs.utility.TimeUtility;
import static io.github.malapert.jwcs.utility.TimeUtility.epochs;


/**
 * Mean place post IAU 1976 system. 
 * 
 * <p>Also a barycentric equatorial coordinate system. 
 * This should be qualified by an Equinox value (only).
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class FK5 implements CoordinateReferenceFrame {    
    /**
     * The name of this reference frame.
     */
    private static final CoordinateReferenceFrame.ReferenceFrame REF_SYSTEM = CoordinateReferenceFrame.ReferenceFrame.FK5;
    
    /**
     * The default value of the equinox sets to J2000.
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
     * Creates a FK5 reference frame with the {@link FK5#DEFAULT_EPOCH} value.
     */
    public FK5(){
        init(DEFAULT_EPOCH);
    }
    
    /**
     * Creates a FK5 reference frame with a equinox value.
     * 
     * <p>The epoch of equinox is transformed in Julian epoch
     * using {@link TimeUtility#epochs(java.lang.String)}. 
     * 
     * @param epoch the epoch of equinox
     */    
    public FK5(final String epoch) {
        init(epoch);
    } 

    /**
     * Initialization.
     * @param epoch the epoch of equinox
     */
    private void init(final String epoch) {          
        this.setEquinox(epoch);
    }

    /**
     * Returns the equinox as a Julian epoch.
     */    
    @Override
    public double getEquinox() {
        return this.equinox;
    }

    /**
     * Returns Double.NaN.
     * 
     * <p>No need to specify an epoch of observation in FK5 reference frame
     */    
    @Override
    public double getEpochObs() {
        return Double.NaN;
    }    

    @Override
    public ReferenceFrame getReferenceFrame() {
        return REF_SYSTEM;
    }

    /**
     * Sets the equinox and transforms it in a Julian epoch.
     * 
     * <p>The epoch of equinox is transformed in Julian epoch
     * using {@link TimeUtility#epochs(java.lang.String) }. 
     * 
     * @param equinox the equinox to set
     */
    @Override
    public final void setEquinox(final String equinox) {
        this.equinox = epochs(equinox)[1];
    }
    
    @Override
    /**
     * Sets the equinox as a Julian epoch.
     */
    public void setEquinox(final double equinox) {
        this.equinox = equinox;
    }      
    
    @Override
    /**
     * Do nothing.
     * 
     * <p>No need to specify the epoch of observation in FK5 reference frame
     */
    public void setEpochObs(final String epochObs) { 
        //Do nothing 
    }

    @Override
    /**
     * Do nothing.
     * 
     * <p>No need to specify the epoch of observation in FK5 reference frame
     */    
    public void setEpochObs(final double epochObs) {
        //Intentionally empty
    }  

    @Override
    public String toString() {
        return "FK5(J"+this.equinox+")";
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 29 * hash + (int) (Double.doubleToLongBits(this.equinox) ^ (Double.doubleToLongBits(this.equinox) >>> 32));
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
        final FK5 other = (FK5) obj;
        return Double.doubleToLongBits(this.equinox) == Double.doubleToLongBits(other.equinox);
    }
        
}
