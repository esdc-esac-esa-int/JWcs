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
 * The old FK4 (barycentric) equatorial system
 * but without the <b>E-terms of aberration</b>.
 * 
 * <p>This coordinate system should also be 
 * qualified by both an Equinox and an Epoch value.
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class FK4NoEterms implements CoordinateReferenceFrame {

    /**
     * The name of this reference frame.
     */
    private static final CoordinateReferenceFrame.ReferenceFrame REF_SYSTEM = CoordinateReferenceFrame.ReferenceFrame.FK4_NO_E;
    
    /**
     * The default value of the epoch sets to B1950.
     */
    private static final String DEFAULT_EPOCH = "B1950";
    
    /**
     * The epoch of the equinox.
     * 
     * <p>An equinox is an astronomical event in which the plane of Earth's 
     * equator passes through the center of the Sun, which occurs twice each 
     * year, around 20 March and 23 September.     
     */
    private double equinox;
    
    /**
     * The epoch of observation.
     */
    private double epochObs;
    
    /**
     * Creates a FK4_NO_E reference frame with default value of the epoch of equinox {@link FK4NoEterms#DEFAULT_EPOCH}.
     */
    public FK4NoEterms() {
        init(DEFAULT_EPOCH, DEFAULT_EPOCH);
    }
    
    /**
     * Creates a FK4_NO_E reference frame with a equinox.
     * @param epoch the epoch of equinox
     */    
    public FK4NoEterms(final String epoch) {       
        init(epoch, null);
    }
    
    /**
     * Creates a FK4NoEterms reference frame with both equinox and epoch of observation.
     * @param epoch the epoch
     * @param epochObs the epoch of observation
     */    
    public FK4NoEterms(final String epoch, final String epochObs) {       
        init(epoch, epochObs);
    }
    
    /**
     * Initialization.
     * @param epoch the epoch of equinox
     * @param epochObs the epoch of observation
     */
    private void init(final String epoch, final String epochObs) {        
        this.setEpochObs(epochObs);        
        this.setEquinox(epoch);
    }
    
    @Override
    public ReferenceFrame getReferenceFrame() {
        return REF_SYSTEM;
    }

    @Override
    /**
     * Returns the Besselian value of the epoch of observation.
     */    
    public double getEpochObs() {
        return this.epochObs;
    }

    @Override
    /**
     * Returns the Besselian value of the equinox.
     */    
    public double getEquinox() {
        return this.equinox;
    }

    /**
     * Sets the equinox.
     * 
     * <p>The epoch of equinox is transformed in Besselian epoch
     * using {@link TimeUtility#epochs(java.lang.String) }. 
     * 
     * @param equinox the equinox to set
     */
    @Override
    public final void setEquinox(final String equinox) {
        this.equinox = epochs(equinox)[0];
    }        

    /**
     * Sets the epoch of observation.
     * 
     * <p>The epoch of observation is transformed in Besselian epoch
     * using {@link TimeUtility#epochs(java.lang.String) }. 
     * 
     * <p>When epochObs is null, the epochObs is set to NaN.
     * 
     * @param epochObs the epochObs to set
     */
    @Override
    public final void setEpochObs(final String epochObs) {
        this.epochObs = (epochObs == null) ? Double.NaN : epochs(epochObs)[0];
    }
    
    /**
     * Sets the Besselian epoch of the equinox.
     * @param equinox the equinox to set
     */
    @Override
    public void setEquinox(final double equinox) {
        this.equinox = equinox;
    }        

    /**
     * Sets the Besselian epoch of the epoch of observation.
     * @param epochObs the epochObs to set
     */
    @Override
    public void setEpochObs(final double epochObs) {
        this.epochObs = epochObs;
    }    
    
    @Override
    public String toString() {
        final String result;
        if (Double.isNaN(this.epochObs)) {
            result = "FK4_NO_E(B"+this.equinox+")";
        } else {
            result = "FK4_NO_E(B"+this.equinox+",B"+this.epochObs+")";
        }
        return result;
    } 

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 97 * hash + (int) (Double.doubleToLongBits(this.equinox) ^ (Double.doubleToLongBits(this.equinox) >>> 32));
        hash = 97 * hash + (int) (Double.doubleToLongBits(this.epochObs) ^ (Double.doubleToLongBits(this.epochObs) >>> 32));
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
        final FK4NoEterms other = (FK4NoEterms) obj;
        if (Double.doubleToLongBits(this.equinox) != Double.doubleToLongBits(other.equinox)) {
            return false;
        }
        return Double.doubleToLongBits(this.epochObs) == Double.doubleToLongBits(other.epochObs);
    }
    
}
