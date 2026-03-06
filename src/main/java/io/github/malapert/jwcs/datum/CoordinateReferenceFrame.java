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

import io.github.malapert.jwcs.proj.exception.JWcsError;
import io.github.malapert.jwcs.utility.TimeUtility;

/**
 * The coordinate reference frame defines how the CRS is related to the origin
 * (position and the date of the origin - equinox, date of observation ). 
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 * @see io.github.malapert.jwcs.crs.AbstractCrs
 */
public interface CoordinateReferenceFrame {

    /**
     * The supported reference frames.
     * 
     * <p>The coordinate reference frame or reference frame defines how the CRS 
     * is related to the origin (position and the date of the origin - equinox, 
     * date of observation).
     */
    enum ReferenceFrame {
        /**
         * The International Celestial Reference System, for optical data 
         * realized through the Hipparcos catalog. 
         * 
         * <p>By definition, ICRS is not an equatorial system, but it is 
         * very close to the FK5 (J2000) system. No Equinox value is required.
         */
        ICRS("ICRS", false, false),
        /**
         * Mean place post IAU 1976 system. 
         * 
         * <p>Also a barycentric equatorial coordinate system. 
         * This should be qualified by an Equinox value (only).
         */
        FK5("FK5", true, false), 
        /**
         * Mean place pre-IAU 1976 system. 
         * 
         * <p>FK4 is the old barycentric (i.e. w.r.t. the common 
         * center of mass) equatorial coordinate 
         * system, which should be qualified by an Equinox value.
         * For accurate work FK4 coordinate systems should also be qualified
         * by an Epoch value. This is the *epoch of observation*.
         * 
         */
        FK4("FK4", true, true),
        /**
         * The old FK4 (barycentric) equatorial system
         * but without the *E-terms of aberration*.
         * 
         * <p>This coordinate system should also be 
         * qualified by both an Equinox and an Epoch value.
         * 
         * <p>In 'Representations of celestial coordinates in FITS' (Calabretta and Greisen) 
         * we read that all reference systems are allowed for both equatorial and 
         * ecliptic coordinates, except FK4-NO-E, which is only allowed for equatorial 
         * coordinates. If FK4-NO-E is given in combination with an ecliptic 
         * crs then silently FK4 is assumed.         
         * 
         * @see <a href="http://www.atnf.csiro.au/people/mcalabre/WCS/ccs.pdf">
         * "Representations of celestial coordinates in FITS, M. R. Calabretta and E. W. Greisen - page 6 "</a> 
         *          
         */
        FK4_NO_E("FK4 NO E-terms", true, true), 
        /**
         * This is an equatorial coordinate system based on 
         * the mean dynamical equator and equinox at epoch J2000.
         * 
         * <p>The dynamical equator and equinox differ slightly 
         * compared to the equator and equinox of FK5 at J2000 and 
         * the ICRS system. This system need not be qualified 
         * by an Equinox value.
         */
        J2000("J2000", false, false);
        
        /**
         * Name of the reference frame.
         */
        private final String name;
        
        /**
         * Needs an equinox value as parameter.
         */
        private final boolean hasEquinox;
        
        /**
         * Needs an epoch value as parameter.
         */
        private final boolean hasEpoch;
        
        /**
         * Constructor.
         * @param name Name of the reference frame
         * @param hasEquinox Can have equinox as parameter
         * @param hasEpoch Can have epoch as parameter
         */
        ReferenceFrame(final String name, final boolean hasEquinox, final boolean hasEpoch) {
            this.name = name;
            this.hasEquinox = hasEquinox;
            this.hasEpoch = hasEpoch;
        }
        
        /**
         * Returns the name of the CRS.
         * @return the name of the CRS
         */
        public String getName() {
            return this.name;
        }
        
        /**
         * Returns True when the reference frame needs the equinox as parameter.
         * @return True when the reference frame needs the equinox as parameter otherwise False
         */
        public boolean hasEquinox() {
            return this.hasEquinox;
        }

        /**
         * Returns True when the reference frame needs the epoch as parameter.
         * @return True when the reference frame needs the epoch as parameter otherwise False
         */        
        public boolean hasEpoch() {
            return this.hasEpoch;
        }
        
        /**
         * Returns the CoordinateReferenceFrame based on its name.
         * @param name name of the reference frame
         * @return the CoordinateReferenceFrame type
         * @throws JWcsError Reference frame not found by searching by its name
         */
        public static ReferenceFrame valueOfByName(final String name) {
            ReferenceFrame result = null;
            final ReferenceFrame[] values = ReferenceFrame.values();
            for (ReferenceFrame value : values) {
                if(value.getName().equals(name)) {
                    result = value;
                    break;
                }
            }
            if (result == null) {
                throw new JWcsError("Reference frame not found by searching by its name "+name);
            } else {
                return result;
            }
        }
        
        /**
         * Returns the list of CoordinateReferenceFrame.
         * @return all names of the CoordinateReferenceFrame
         */
        public static String[] getRefenceFrameNametoArray() {            
            final ReferenceFrame[] values = ReferenceFrame.values();
            final String[] result = new String[values.length];
            int index = 0;
            for (ReferenceFrame value : values) {
                result[index] = value.getName();
                index++;
            }
            return result;
        }                
    }
   
    /**
     * Returns the reference frame that is used.
     * @return the reference frame that is used
     */
    CoordinateReferenceFrame.ReferenceFrame getReferenceFrame();
    
    /**
     * Returns the epoch of observation as a Besselian or a Julian value 
     * according to the reference frame.
     * @return Double.NaN when epoch of observation is not required otherwise
     * the epoch of observation
     * @see TimeUtility#convertEpochBessel2JD
     * @see TimeUtility#convertEpochJulian2JD
     */
    double getEpochObs();
    
    /**
     * Returns the equinox as a Besselian or a Julian value according to the 
     * reference frame.
     * @return the equinox
     * @see TimeUtility#convertEpochBessel2JD
     * @see TimeUtility#convertEpochJulian2JD
     */
    double getEquinox();
    
    /**
     * Sets the epoch of observation.
     * @param epochObs the epoch of observation
     * @see TimeUtility#epochs
     */
    void setEpochObs(final String epochObs);
    
    /**
     * Sets the equinox.
     * @param equinox the equinox
     * @see TimeUtility#epochs
     */
    void setEquinox(final String equinox);
    
    /**
     * Sets the epoch of observation as a Julian or Besselian epoch according to
     * the reference frame.
     * @param epochObs the epoch of observation
     * @see TimeUtility#convertEpochBessel2JD
     * @see TimeUtility#convertEpochJulian2JD     
     */
    void setEpochObs(final double epochObs);
    
    /**
     * Sets the equinox as a Julian or Besselian epoch according to
     * the reference frame.
     * @param equinox the equinox
     * @see TimeUtility#convertEpochBessel2JD
     * @see TimeUtility#convertEpochJulian2JD
     */
    void setEquinox(final double equinox);    
       
}
