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
 * This is an equatorial coordinate reference frame based on 
 * the mean dynamical equator and equinox at epoch J2000.
 * 
 * <p>The dynamical equator and equinox differ slightly 
 * compared to the equator and equinox of FK5 at J2000 and 
 * the ICRS system. This system need not be qualified 
 * by an Equinox value.
 *         
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class J2000 implements CoordinateReferenceFrame {

    /**
     * The name of this reference frame.
     */
    private static final CoordinateReferenceFrame.ReferenceFrame REF_SYSTEM = CoordinateReferenceFrame.ReferenceFrame.J2000;

    /**
     * The default value of the equinox sets to J2000.
     */
    private static final String DEFAULT_EPOCH = "J2000";

    /**
     * The epoch of the equinox.
     */
    private double equinox;

    /**
     * Creates J2000 frame based on {@link J2000#DEFAULT_EPOCH}.
     */
    public J2000() {
        init(DEFAULT_EPOCH);
    }
    
    /**
     * Initialization.
     *
     * @param epoch the epoch of equinox
     */
    private void init(final String epoch) {
        this.equinox = epochs(epoch)[1];
    }

    @Override
    public ReferenceFrame getReferenceFrame() {
        return REF_SYSTEM;
    }

    /**
     * Returns NaN.
     * 
     * <p>No need to specify an epoch of observation in J2000 reference frame
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
     * Sets the equinox and transforms it as an Julian epoch.
     * @see TimeUtility#epochs(java.lang.String) 
     */ 
    @Override   
    public void setEquinox(final String equinox) { 
        this.equinox = epochs(equinox)[1];
    }    
    
    /**
     * Sets the equinox as Julian epoch.
     */    
    @Override
    public void setEquinox(final double equinox) {
        this.equinox = equinox;
    }      
    
    /**
     * Do nothing.
     * 
     * <p>No need to specify the epoch of observation
     */    
    @Override
    public void setEpochObs(final String epochObs) { 
        //Do nothing
    }

    /**
     * Do nothing.
     * 
     * <p>No need to specify the epoch of observation
     */     
    @Override   
    public void setEpochObs(final double epochObs) {
        //Intentionally empty because the epoch is invariant in J2000 CRS.
    }     

    @Override
    public String toString() {
        return "J2000";
    }        

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 97 * hash + (int) (Double.doubleToLongBits(this.equinox) ^ (Double.doubleToLongBits(this.equinox) >>> 32));
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
        final J2000 other = (J2000) obj;
        return Double.doubleToLongBits(this.equinox) == Double.doubleToLongBits(other.equinox);
    }

}
