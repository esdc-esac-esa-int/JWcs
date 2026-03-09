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

import static io.github.malapert.jwcs.utility.NumericalUtility.createRealMatrix;
import io.github.malapert.jwcs.utility.TimeUtility;
import static io.github.malapert.jwcs.utility.TimeUtility.epochs;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

/**
 * Mean place pre-IAU 1976 system.
 * 
 * <p>FK4 is the old barycentric (i.e. w.r.t. the common center of mass) equatorial
 * coordinate system, which should be qualified by an Equinox value.
 * For accurate work FK4 coordinate systems should also be qualified
 * by an Epoch value. This is the <b>epoch of observation</b>.
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class FK4 implements CoordinateReferenceFrame {
    /**
     * The name of this reference frame.
     */
    private static final CoordinateReferenceFrame.ReferenceFrame REF_SYSTEM = CoordinateReferenceFrame.ReferenceFrame.FK4;
    
    /**
     * The default value of the epoch set to B1950.
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
     * 
     * <p>This is the time when the mean place was correct.
     */    
    private double epochObs;       
    
    /**
     * Creates a FK4 reference frame with default value of epoch {@link FK4#DEFAULT_EPOCH}.
     */
    public FK4() {
        init(DEFAULT_EPOCH, DEFAULT_EPOCH);
    }
    
    /**
     * Creates a FK4 reference frame with an epoch of equinox.
     * @param epoch the epoch of equinox
     */
    public FK4(final String epoch) {       
        init(epoch, epoch);
    }
    
    /**
     * Creates a FK4 reference frame with both equinox and epoch of observation.
     * 
     * <p>The epoch of equinox and epoch of observation are transformed in 
     * Besselian epochs using {@link TimeUtility#epochs(java.lang.String)}.
     * 
     * @param epoch the epoch of equinox
     * @param epochObs the epoch of observation
     * @see TimeUtility#epochs(java.lang.String) 
     */
    public FK4(final String epoch, final String epochObs) {       
        init(epoch, epochObs);
    }
    
    /**
     * initialization.
     * 
     * <p>The epoch of equinox and epoch of observation are transformed in 
     * Besselian epochs using {@link TimeUtility#epochs(java.lang.String)}. 
     * 
     * @param epoch the epoch of equinox
     * @param epochObs the epoch of observation
     * @see TimeUtility#epochs(java.lang.String)      
     */
    private void init(final String epoch, final String epochObs) {        
        this.setEpochObs(epochObs);
        this.setEquinox(epoch);
    }
    
    /**
     * Compute the E-terms (elliptic terms of aberration) for a given epoch.
     * 
     * <p>Reference:<br>
     * ----------<br>
     * Seidelman, P.K.,  1992.  Explanatory Supplement to the Astronomical
     * Almanac.  University Science Books, Mill Valley
     *
     * <p>Notes:<br>     
     * -------<br>
     * The method is described on page 170/171 of the ES.
     * One needs to process the e-terms for the appropriate
     * epoch This routine returns the e-term vector for arbitrary epoch.
     * 
     * @param epoch A Besselian epoch
     * @return A tuple containing the e-terms vector (DeltaD,DeltaC,DeltaC.tan(e0))
     */
    public static RealMatrix getEterms(final double epoch) {
        //Julian centuries since B1950
        final double T = (epoch-1950.0d)*1.00002135903d/100.0d;
        //Eccentricity of the Earth's orbit
        final double ec = 0.01673011d-(0.00004193d+0.000000126d*T)*T;
        //Mean obliquity of the ecliptic. Method is different compared to 
        //functions for the obliquity defined earlier. This function depends
        //on time wrt. epoch 1950 not epoch 2000.
        double ob = 84404.836d-(46.8495d+(0.00319d+0.00181d*T)*T)*T;
        ob = FastMath.toRadians(ob/3600.0d);
        //Mean longitude of perihelion of the solar orbit
        double p = 1015489.951d+(6190.67d+(1.65d+0.012d*T)*T)*T;
        p = FastMath.toRadians(p/3600.0d);
        //Calculate the E-terms vector
        final double ek = ec*FastMath.toRadians(20.49522d/3600.0d); // 20.49552 is constant of aberration at J2000
        final double cp = FastMath.cos(p);
        //       -DeltaD        DeltaC            DeltaC.tan(e0)
        final double[][] array = {
            {ek*FastMath.sin(p), -ek*cp*FastMath.cos(ob), -ek*cp*FastMath.sin(ob)}
        };
        return createRealMatrix(array);       
    } 

    /**
     * Returns the Besselian value of the epoch of observation.
     */    
    @Override
    public double getEpochObs() {
        return this.epochObs;
    }

    /**
     * Returns the Besselian epoch of the equinox.
     */    
    @Override
    public double getEquinox() {
        return this.equinox;
    }   

    @Override
    public ReferenceFrame getReferenceFrame() {
        return REF_SYSTEM;
    }

    /**
     * Sets the equinox and transforms it in a Besselian epoch.
     * 
     * <p>The epoch of equinox is transformed in Besselian epochs 
     * using {@link TimeUtility#epochs(java.lang.String)}. 
     * 
     * @param equinox the equinox to set
     * @see TimeUtility#epochs(java.lang.String) 
     */
    @Override
    public final void setEquinox(final String equinox) {
        this.equinox = epochs(equinox)[0];
    }
    
    /**
     * Sets the Besselian epoch of the equinox.
     * @param equinox the equinox to set
     */
    @Override
    public final void setEquinox(final double equinox) {
        this.equinox = equinox;
    }    

    /**
     * Sets the epoch and transforms it in a Besselian epoch.
     * @param epochObs the epochObs to set
     */
    @Override
    public final void setEpochObs(final String epochObs) {
        this.epochObs = (epochObs == null)? Double.NaN : epochs(epochObs)[0];
    }
    
    /**
     * Sets the Besselian epoch of observation.
     *    
     * @param epochObs the epochObs to set
     */
    @Override
    public final void setEpochObs(final double epochObs) {
        this.epochObs = epochObs;
    }    

    @Override
    public String toString() {
        return Double.isNaN(this.epochObs) ? "FK4(B"+this.equinox+")" : "FK4(B"+this.equinox+", B"+this.epochObs+")";
    }    

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 17 * hash + (int) (Double.doubleToLongBits(this.equinox) ^ (Double.doubleToLongBits(this.equinox) >>> 32));
        hash = 17 * hash + (int) (Double.doubleToLongBits(this.epochObs) ^ (Double.doubleToLongBits(this.epochObs) >>> 32));
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
        final FK4 other = (FK4) obj;
        if (Double.doubleToLongBits(this.equinox) != Double.doubleToLongBits(other.equinox)) {
            return false;
        }
        return Double.doubleToLongBits(this.epochObs) == Double.doubleToLongBits(other.epochObs);
    }
    
}
