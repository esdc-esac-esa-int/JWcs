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
package io.github.malapert.jwcs.position;

import io.github.malapert.jwcs.crs.AbstractCrs;
import io.github.malapert.jwcs.utility.DMS;
import io.github.malapert.jwcs.utility.HMS;
import static io.github.malapert.jwcs.utility.NumericalUtility.aacos;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.util.FastMath;

/**
 * Represents a position in the sky.
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class SkyPosition {
    
    /**
     * Logger.
     */
    private static final Logger LOG = Logger.getLogger(SkyPosition.class.getName());
    
    
    /**
     * Factor to convert degrees to hours.
     */
    public static final double TO_HOURS = 24d/360d;
    
    /**
     * Longitude in decimal degrees of the position.
     * @see #getLongitude() 
     * @see #getLongitudeAsSexagesimal() 
     */
    private double longitude;
    /**
     * Latitude in decimal degrees of the position.
     * @see #getLatitude() 
     * @see #getLatitudeAsSexagesimal()  
     */    
    private double latitude;
    /**
     * AbstractCrs of the position.
     * @see #getCrs() 
     */
    private AbstractCrs crs;

    /**
     * Creates a position in the sky.
     * 
     * <p>The position can be created in different kind of CRS.
     * 
     * @param longitude longitude in decimal degrees of the position
     * @param latitude latitude in decimal degrees of the position
     * @param crs CRS in which the position is stored
     */
    public SkyPosition(final double longitude, final double latitude, final AbstractCrs crs) {
        this.longitude = longitude%360;
        this.latitude = latitude;
        this.crs = crs;
    } 

    /**
     * Returns the longitude in decimal degrees.
     * @return the longitude
     */
    public double getLongitude() {
        return longitude;
    }
    
    /**
     * Returns the longitude in sexagesimal HH:MM:SS.SS.
     * @return the longitude in sexagesimal HH:MM:SS.SS
     */
    public String getLongitudeAsSexagesimal() {
        final HMS hms = new HMS(getLongitude()*TO_HOURS);
        return hms.toString(true);
    }

    /**
     * Sets the longitude in decimal degrees.
     * @param longitude the longitude to set
     */
    public void setLongitude(final double longitude) {
        this.longitude = longitude;
    }

    /**
     * Returns the latitude in decimal degrees.
     * @return the latitude
     */
    public double getLatitude() {
        return latitude;
    }

    /**
     * Sets the latitude in decimal degrees.
     * @param latitude the latitude to set
     */
    public void setLatitude(final double latitude) {
        this.latitude = latitude;
    }
    
    /**
     * Returns the latitude in sexagesimal DD:MM:SS.SS.
     * @return the latitude
     */
    public String getLatitudeAsSexagesimal() {
        final DMS dms = new DMS(getLatitude());
        return dms.toString(true);
    }   
    
    /**
     * Returns the skyposition as a double array.
     * @return the sky position.
     */
    public double[] getDoubleArray() {
        return new double[]{getLongitude(), getLatitude()};
    }

    /**
     * Returns the coordinate reference system.
     * @return the crs
     */
    public AbstractCrs getCrs() {
        return crs;
    }

    /**
     * Sets the coordinate reference system.
     * @param crs the crs to set
     */
    public void setCrs(final AbstractCrs crs) {
        this.crs = crs;
    }

    /**
     * Returns the cartesian coordinates in the current reference frame.
     * @return the cartesian coordinates as an array
     */
    public double[] getCartesian() {
        return AbstractCrs.longlat2xyz(this.longitude, this.latitude).getColumn(0);
    }
    
    /**
     * Computes the angular separation between two positions in different
     * coordinate reference systems.
     *
     * @param pos1 sky position in a coordinate Reference System
     * @param pos2 sky position in a coordinate Reference System
     * @return angular separation in decimal degrees.
     */
    public static double separation(final SkyPosition pos1, final SkyPosition pos2) {
        final AbstractCrs crs = pos1.getCrs();
        final SkyPosition pos1InRefFramePos2 = crs.convertTo(pos2.getCrs(), pos1.getLongitude(), pos1.getLatitude());
        final double[] pos1XYZ = pos1InRefFramePos2.getCartesian();
        final double[] pos2XYZ = pos2.getCartesian();
        final double normPos1 = FastMath.sqrt(pos1XYZ[0] * pos1XYZ[0] + pos1XYZ[1] * pos1XYZ[1] + pos1XYZ[2] * pos1XYZ[2]);
        final double normPos2 = FastMath.sqrt(pos2XYZ[0] * pos2XYZ[0] + pos2XYZ[1] * pos2XYZ[1] + pos2XYZ[2] * pos2XYZ[2]);
        final double separation = aacos((pos1XYZ[0] * pos2XYZ[0] + pos1XYZ[1] * pos2XYZ[1] + pos1XYZ[2] * pos2XYZ[2]) / (normPos1 * normPos2));
        LOG.log(Level.INFO, "seratation({0},{1}) =  {2}", new Object[]{pos1, pos2, FastMath.toDegrees(separation)});
        return FastMath.toDegrees(separation);
    }     
      
    @Override
    public String toString() {        
        return "SkyPosition : (" + getLongitude() + "," + getLatitude() + ") or (" + getLongitudeAsSexagesimal() + "," + getLatitudeAsSexagesimal() + ") in " + getCrs();
    }  

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 53 * hash + (int) (Double.doubleToLongBits(this.longitude) ^ (Double.doubleToLongBits(this.longitude) >>> 32));
        hash = 53 * hash + (int) (Double.doubleToLongBits(this.latitude) ^ (Double.doubleToLongBits(this.latitude) >>> 32));
        hash = 53 * hash + Objects.hashCode(this.crs);
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
        final SkyPosition other = (SkyPosition) obj;
        if (Double.doubleToLongBits(this.longitude) != Double.doubleToLongBits(other.longitude)) {
            return false;
        }
        if (Double.doubleToLongBits(this.latitude) != Double.doubleToLongBits(other.latitude)) {
            return false;
        }
        return Objects.equals(this.crs, other.crs);
    }
    
    
}
