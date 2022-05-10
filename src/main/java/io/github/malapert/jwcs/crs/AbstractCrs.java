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

import io.github.malapert.jwcs.position.SkyPosition;
import io.github.malapert.jwcs.datum.FK4;
import io.github.malapert.jwcs.datum.CoordinateReferenceFrame;
import io.github.malapert.jwcs.proj.exception.JWcsError;
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.aasin;
import static io.github.malapert.jwcs.utility.NumericalUtility.createRealIdentityMatrix;
import static io.github.malapert.jwcs.utility.NumericalUtility.createRealMatrix;
import static io.github.malapert.jwcs.utility.NumericalUtility.inverse;
import static io.github.malapert.jwcs.utility.NumericalUtility.rotX;
import static io.github.malapert.jwcs.utility.NumericalUtility.rotY;
import static io.github.malapert.jwcs.utility.NumericalUtility.rotZ;
import static io.github.malapert.jwcs.utility.TimeUtility.convertEpochBessel2JD;
import static io.github.malapert.jwcs.utility.TimeUtility.convertEpochJulian2JD;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import io.github.malapert.jwcs.utility.TimeUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.aatan2;
import static io.github.malapert.jwcs.utility.NumericalUtility.equal;
import static io.github.malapert.jwcs.utility.NumericalUtility.isInInterval;

/**
 * A Coordinate Reference System (crs) contains two different elements : 
 * the <b>coordinate reference frame</b> {@link AbstractCrs#getCoordinateReferenceFrame()}
 * and the <b>coordinate system</b> {@link AbstractCrs#getCoordinateReferenceSystem()} .
 *
 * <p>The coordinate reference frame defines how the CRS is related to the origin
 * (position and the date of the origin - equinox {@link CoordinateReferenceFrame#getEquinox()} , 
 * date of observation {@link CoordinateReferenceFrame#getEpochObs()} ) and 
 * the coordinate system describes how the coordinates expressed in the 
 * coordinate reference frame (e.g. as cartesian coordinates, spherical 
 * coordinates or coordinates of a map projection).
 * 
 * <p>An equinox is an astronomical event in which the plane of Earth's equator 
 * passes through the center of the Sun.
 * 
 * <p>An epoch of observation is the moment in time when the coordinates are known 
 * to be correct. Often, this will be the date of observation, and is important 
 * in cases where coordinates systems move with respect to each other over the 
 * course of time
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public abstract class AbstractCrs {

    /**
     * Logger.
     */
    private static final Logger LOG = Logger.getLogger(AbstractCrs.class.getName());

    /**
     * List of supported CoordinateReferenceSystem. 
     */
    public enum CoordinateReferenceSystem {
        /**
         * Galactic coordinates (lII, bII).
         */
        GALACTIC("Galactic", false),
        /**
         * Equatorial coordinates (\u03B1, \u03B4).
         */
        EQUATORIAL("Equatorial", true),
        /**
         * De Vaucouleurs Supergalactic coordinates (sgl, sgb).
         */
        SUPER_GALACTIC("Super Galactic", false),
        /**
         * Ecliptic coordinates (\u03BB, \u03B2) referred to the ecliptic and
         * mean equinox.
         */
        ECLIPTIC("Ecliptic", true);
        
        /**
         * Coordinate system name.
         */
        private final String name;
        /**
         * Indicates if the coordinate system has a coordinate reference frame.
         */
        private final boolean hasCoordinateReferenceFrame;
        
        /**
         * Constructor.
         * @param name coordinate system name
         * @param hasCoordinateReferenceFrame Indicates if the coordinate system has a 
         * configurable coordinate reference frame.
         */
        CoordinateReferenceSystem(final String name, final boolean hasCoordinateReferenceFrame) {
            this.name = name;
            this.hasCoordinateReferenceFrame = hasCoordinateReferenceFrame;
        }
        
        /**
         * Returns the name of the coordinate system.
         * @return the coordinate system
         */
        public final String getName() {
            return this.name;
        }
        
        /**
         * Tests if the the coordinate system has a coordinate reference frame.
         * @return True when the coordinate system has a coordinate reference 
         * frame otherwise false
         */
        public final boolean hasCoordinateReferenceFrame() {
            return this.hasCoordinateReferenceFrame;
        }
      
        /**
         * Returns the CoordinateReferenceSystem object based on its name.
         * @param name name of the coordinate system
         * @return the CoodinateSystem object
         * @exception JWcsError Coordinate system not found
         */
        public static CoordinateReferenceSystem valueOfByName(final String name) {
            CoordinateReferenceSystem result = null;
            final CoordinateReferenceSystem[] values = CoordinateReferenceSystem.values();
            for (final CoordinateReferenceSystem value : values) {
                if(value.getName().equals(name)) {
                    result = value;
                    break;
                }
            }
            if (result == null) {
                throw new JWcsError("Coordinate system not found by searching by its name "+name);
            } else {
                return result;
            }
        }
        
        /**
         * Returns an array of the coordinate system name.
         * @return An array of names of all coordinate system
         */
        public static String[] getCoordinateSystemArray() {            
            final CoordinateReferenceSystem[] values = CoordinateReferenceSystem.values();
            String[] result = new String[values.length];
            int index = 0;
            for (final CoordinateReferenceSystem value : values) {
                result[index] = value.getName();
                index++;
            }
            return result;
        }
    }

    /**
     * Computes the rotation matrix to apply from a reference frame to another one.
     *
     * @param crs the output coordinate Reference System
     * @return the rotation matrix in the output coordinate Reference System
     * @throws JWcsError Unknown output crs
     */        
    protected abstract RealMatrix getRotationMatrix(final AbstractCrs crs) throws JWcsError;

    /**
     * Returns the coordinate system of the CRS.     
     * @return the type of the coordinate system
     */
    public abstract CoordinateReferenceSystem getCoordinateReferenceSystem();    
    
    /**
     * Returns the coordinate reference frame of the CRS.
     * @return the coordinate reference frame of the CRS
     */
    public abstract CoordinateReferenceFrame getCoordinateReferenceFrame();    

    /**
     * Sets the coordinate reference frame of the CRS.
     * @param refFrame the coordinate reference frame
     */
    public abstract void setCoordinateReferenceFrame(final CoordinateReferenceFrame refFrame);      
    
    /**
     * Returns the elliptical terms of aberration matrix for the input 
     * coordinate Reference System.
     * 
     * <p>Stellar aberration is caused by the motion of the earth in its orbit. 
     * This motion is represented by a circular velocity component and a 
     * component perpendicular to the major axis caused by the fact that the 
     * orbit is elliptical. This velocity component is responsible for 
     * elliptical terms of aberration (E-terms) which are less than 
     * 0.35 arcseconds (maximum is equal to the constant of aberration times 
     * the eccentricity of the earths orbit = 20”.496 x 0.01673 ~= 343 mas). 
     * The terms are independent of the position of the earth and depend only 
     * on the position of the object in the sky.
     * 
     * <p>The elliptical terms of aberration are only applied for FK4 coordinate
     * reference frame.
     *
     * @param crs the coordinate reference system
     * @return the elliptical terms of aberration matrix
     * @see FK4#getEterms(double) 
     */
    protected static final RealMatrix getEterms(final AbstractCrs crs) {
        RealMatrix eterms = null;
        final CoordinateReferenceFrame refFrame = crs.getCoordinateReferenceFrame();        
        if (refFrame!=null && CoordinateReferenceFrame.ReferenceFrame.FK4.equals(refFrame.getReferenceFrame())) {
            final double equinox = refFrame.getEquinox();
            eterms = FK4.getEterms(equinox);
            LOG.log(Level.FINER, "getEterms {0}(FK4) from {1} : {2}", new Object[]{crs.getCoordinateReferenceSystem(),equinox,eterms});
        }
        return eterms;        
    }
  
    /**
     * Checks the coordinates.
     * 
     * <p>The longitude must be included in [0,360] while the latitude must be 
     * included in [-90,90].
     * 
     * @param longitude longitude in decimal degrees
     * @param latitude  latitude in decimal degrees
     * @throws JWcsError latitude or longitude is out of range
     */
    private void checkCoordinates(final double longitude, final double latitude) throws JWcsError {
        final boolean isLongInterval = isInInterval(longitude, 0, 360);
        final boolean isLatInterval = isInInterval(latitude, -90, 90);        
        if(!isLongInterval && !isLatInterval) {
            throw new JWcsError("longitude must be in [0,360] and latitude in [-90,90]");
        } else if (!isLongInterval) {
            throw new JWcsError("longitude must be in [0,360]");    
        } else if (!isLatInterval) {
            throw new JWcsError("latitude must be in [-90,90]");    
        }
    }

    /**
     * Converts the (longitude, latitude) coordinates from the current 
     * coordinate reference system into the target coordinate reference System.     
     * 
     * <p>The algorithm to convert coordinates from a CRS to another is the following:
     * <ul>
     * <li>Transforms coordinates to rectangular coordinates</li>
     * <li>Removes E-terms of source CRS from rectangular coordinates when necessary</li>
     * <li>Applies the rotation matrix from the target CRS to rectangular coordinates</li>
     * <li>Adds E-terms of target CRS from the previous step when necessary</li>
     * <li>Converts to longitude/latitude</li>
     * </ul>
     *
     * @param crs the target coordinate Reference System
     * @param longitude longitude in decimal degrees
     * @param latitude latitude in decimal degrees
     * @return the position on the sky in the target coordinate reference system
     * @see AbstractCrs#longlat2xyz
     * @see AbstractCrs#getRotationMatrix
     * @see AbstractCrs#getEterms
     * @see AbstractCrs#xyz2longlat
     */
    public final SkyPosition convertTo(final AbstractCrs crs, final double longitude, final double latitude) {
        checkCoordinates(longitude, latitude);
        RealMatrix xyz = longlat2xyz(longitude, latitude);
        LOG.log(Level.FINER, "convert sky ({0},{1}) to xyz : {2}", new Object[]{longitude, latitude, xyz});
        final RealMatrix rotation = getRotationMatrix(crs);
        LOG.log(Level.FINER, "Rotation matrix from {0} to {1} : {2}", new Object[]{this.getCoordinateReferenceSystem(),crs.getCoordinateReferenceSystem(),rotation});        
        final RealMatrix etermsIn = AbstractCrs.getEterms(this);
        LOG.log(Level.FINER, "EtermsIn : {0}", new Object[]{etermsIn});        
        final RealMatrix etermsOut = AbstractCrs.getEterms(crs);
        LOG.log(Level.FINER, "EtermsOut from {0} : {1}", new Object[]{crs.getCoordinateReferenceSystem(), etermsOut});        
        if (etermsIn != null) {
            xyz = removeEterms(xyz, etermsIn);
            LOG.log(Level.FINER, "Remove EtermsIn from xyz : {0}", xyz);            
        }      
        xyz = rotation.multiply(xyz);
        if (etermsOut != null) {
            xyz = addEterms(xyz, etermsOut);
            LOG.log(Level.FINER, "Add EtermsOut to xyz : {0}", xyz);            
        }          
        LOG.log(Level.FINER, "Rotate xyz : {0}", xyz);                    
        final double[] position = xyz2longlat(xyz);
        LOG.log(Level.FINER, "Transforms xyz -> ra,dec : {0},{1}", new Object[]{position[0],position[1]});        
        LOG.log(Level.INFO, "convert ({0},{1}) from {2} to {3} --> ({4},{5})", new Object[]{longitude, latitude, this, crs, position[0], position[1]});
        return new SkyPosition(position[0], position[1], crs);
    }

    /**
     * Converts an array of (longitude1, latitude2, longitude2, latitude2, ...)
     * coordinates from the current coordinate system into the target 
     * coordinate reference system.
     *
     * @param crs the target coordinate reference system
     * @param coordinates an array of (longitude1, latitude2, longitude2,
     * latitude2, ...) in degrees
     * @return an array of SkyPosition
     * @throws JWcsError coordinates should be an array containing a set of [longitude, latitude]
     * @see #convertTo 
     */
    public final SkyPosition[] convertTo(final AbstractCrs crs, final double[] coordinates) throws JWcsError {
        
        final int numberElts = coordinates.length;
        final int numberOfCoordinatesPerPoint = 3;
        if (numberElts % 2 != 0) {
            throw new JWcsError("coordinates should be an array containing a set of [longitude, latitude]");
        }
        final SkyPosition[] skyPositionArray = new SkyPosition[(int) (numberElts * 0.5) * numberOfCoordinatesPerPoint];

        final RealMatrix rotation = getRotationMatrix(crs);
        LOG.log(Level.FINER, "Rotation matrix from {0} to {1} : {2}", new Object[]{this.getCoordinateReferenceSystem(),crs.getCoordinateReferenceSystem(),rotation});
        final RealMatrix etermsIn = AbstractCrs.getEterms(this);
        LOG.log(Level.FINER, "EtermsIn : {0}", etermsIn);
        final RealMatrix etermsOut = AbstractCrs.getEterms(crs);
        LOG.log(Level.FINER, "EtermsOut from {0} : {1}", new Object[]{crs.getCoordinateReferenceSystem(), etermsOut});

        int indice = 0;
        for (int i = 0; i < numberElts; i = i + 2) {
            checkCoordinates(coordinates[i], coordinates[i + 1]);
            RealMatrix xyz = longlat2xyz(coordinates[i], coordinates[i + 1]);
            LOG.log(Level.FINER, "xyz : {0}", xyz);
            if (etermsIn != null) {
                xyz = removeEterms(xyz, etermsIn);
                LOG.log(Level.FINER, "Remove EtermsIn from xyz : {0}", xyz);
            }
            xyz = rotation.multiply(xyz);
            LOG.log(Level.FINER, "Rotate xyz : {0}", xyz);            
            if (etermsOut != null) {
                xyz = addEterms(xyz, etermsOut);
                LOG.log(Level.FINER, "Add EtermsOut to xyz : {0}", xyz);
            }
            final double[] position = xyz2longlat(xyz);
            LOG.log(Level.FINER, "Transforms xyz -> ra,dec : {0}", position);
            skyPositionArray[indice] = new SkyPosition(position[0], position[1], crs);
            indice++;
        }
        LOG.log(Level.INFO, "convert {0} from {1} to {2} --> {3}", new Object[]{Arrays.toString(coordinates), this.getCoordinateReferenceSystem(), crs.getCoordinateReferenceSystem(), Arrays.toString(skyPositionArray)});
        return skyPositionArray;
    }
    
    /**
     * Converts a sky position in a CRS to a target CRS.
     * @param targetCrs target CRS
     * @param position position in a CRS
     * @return a SkyPosition to a target CRS
     */
    public static final SkyPosition convertTo(final AbstractCrs targetCrs, final SkyPosition position) {
        final double longitude = position.getLongitude();
        final double latitude = position.getLatitude();
        final AbstractCrs sourceCrs = position.getCrs();
        return sourceCrs.convertTo(targetCrs, longitude, latitude);
    }    
    
    /**
     * Converts an array of sky positions to a target CRS.
     * @param targetCrs the target CRS
     * @param positions the difference position in different CRS
     * @return an array of sky positions to a target CRS
     */
    public static final SkyPosition[] convertTo(final AbstractCrs targetCrs, final SkyPosition[] positions) {
        final SkyPosition[] targetPositions = new SkyPosition[positions.length];
        int i =0;
        for (final SkyPosition position:positions) {
            targetPositions[i] = convertTo(targetCrs, position);
            i++;
        }
        return targetPositions;
    }
       
    /**
     * It handles precession and the transformation between equatorial
     * reference frames.
     *
     * <p>Notes: Return matrix to transform equatorial coordinates from
     * <code>epoch1</code> to <code>epoch2</code> in either reference frame 
     * FK4 or FK5. Or transforms from epoch, FK4 or FK5 to ICRS or J2000 
     * vice versa. Note that each transformation between FK4 and one of the 
     * other reference systems involves a conversion to FK5 and therefore the 
     * epoch of observation will be involved. 
     *
     * @param epoch1 epoch belonging to coordinate reference frame s1 : either 
     * Besselian or Julian.
     * @param epoch2 epoch belonging to coordinate reference frame s2 : either
     * Besselian or Julian.
     * @param s1 input coordinate reference frame
     * @param s2 output coordinate reference frame
     * @param epobs Epoch of observation. Only valid for conversions between FK4
     * and FK5.
     * @return Rotation matrix to transform a position in one of the coordinate 
     * reference frame <code>s1</code> with <code>epoch1</code> to an equatorial
     * system with equator and equinox at <code>epoch2</code> in reference system <code>S2</code>.
     * @throws JWcsError <code>s1</code> to <code>s2</code> conversion is not supported
     * @see #convertICRSTo 
     * @see #convertFK5To 
     * @see #convertFK4To 
     * @see #convertJ2000To 
     */
    protected static RealMatrix convertMatrixEpoch12Epoch2(final double epoch1, final double epoch2, final CoordinateReferenceFrame.ReferenceFrame s1, final CoordinateReferenceFrame.ReferenceFrame s2, final double epobs) {
       final RealMatrix result;
       switch(s1) {
           case ICRS:
               result = convertICRSTo(s2, epoch2, epobs);
               break;
           case FK5:
               result = convertFK5To(s2, epoch1, epoch2, epobs);
               break;
           case FK4:
           case FK4_NO_E:
               result = convertFK4To(s2, epoch1, epoch2, epobs);
               break;
           case J2000:
               result = convertJ2000To(s2, epoch1, epoch2, epobs);
               break;
           default:
                throw new JWcsError(s1.getName()+" conversion to "+s2.getName()+" is not supported");               
       }
       return result;
    }  

    /**
     * Converts J2000 to another coordinate reference frame.     
     *
     * <p>The matrix to compute J2000 to <code>s2</code> is computed in this way:
     * <ul>
     * <li>for ICRS :  m2.multiply(m1) with m2 = convertICRS2J2000Matrix().transpose()
     * and m1 = convertIAU2006MatrixEpoch12Epoch2(epoch1, 2000.0)</li>
     * <li>for FK5 : m4.multiply(m3).multiply(m2).multiply(m1) with
     * m1 = convertIAU2006MatrixEpoch12Epoch2(epoch1, 2000.0), m2 = convertICRS2J2000Matrix().transpose(),
     * m3 = convertICRS2FK5Matrix() and m4 = convertJulianMatrixEpoch12Epoch2(2000.0, epoch2)</li>
     * <li>for FK4 and FK4-NO-E : m5.multiply(m4).multiply(m3).multiply(m2).multiply(m1)
     * with m1 = convertIAU2006MatrixEpoch12Epoch2(epoch1, 2000.0), m2 = convertICRS2J2000Matrix().transpose(),
     * m3 = convertICRS2FK5Matrix(), m4 = convertFK52FK4Matrix(epobs), 
     * m5 = convertBesselianMatrixEpoch12Epoch2(1950, epoch2)</li>
     * <li>for J2000 : convertIAU2006MatrixEpoch12Epoch2(epoch1, epoch2)</li>
     * </ul>     
     * 
     * @param s2 Output coordinate reference frame
     * @param epoch1 Epoch belonging to coordinate reference frame s1 : either 
     * Besselian or Julian.
     * @param epoch2 Epoch belonging to coordinate reference frame s2 : either 
     * Besselian or Julian.
     * @param epobs Epoch of observation. Only valid for conversions between FK4
     * and FK5.
     * @return Rotation matrix to transform a position in J2000 with <code>epoch1</code>
     * to an equatorial crs with equator and equinox at <code>epoch2</code> in 
     * coordinate reference frame <code>s2</code>.
     * @throws JWcsError J2000 to <code>s2</code> conversion is not supported
     * @see #convertIAU2006MatrixEpoch12Epoch2
     * @see #convertICRS2J2000Matrix
     * @see #convertICRS2FK5Matrix 
     * @see #convertJulianMatrixEpoch12Epoch2
     * @see #convertFK52FK4Matrix
     * @see #convertBesselianMatrixEpoch12Epoch2
     */
    private static RealMatrix convertJ2000To(final CoordinateReferenceFrame.ReferenceFrame s2, final double epoch1, final double epoch2, final double epobs) {
        final RealMatrix result;       
        final RealMatrix m1;
        final RealMatrix m2;
        final RealMatrix m3;        
        final RealMatrix m4;         
        final RealMatrix m5;        
        switch(s2) {
            case ICRS:
                m1 = convertIAU2006MatrixEpoch12Epoch2(epoch1, 2000.0d);
                m2 = convertICRS2J2000Matrix().transpose();
                result = m2.multiply(m1);                              
                break;
            case FK5:
                m1 = convertIAU2006MatrixEpoch12Epoch2(epoch1, 2000.0d);
                m2 = convertICRS2J2000Matrix().transpose();
                m3 = convertICRS2FK5Matrix();
                m4 = convertJulianMatrixEpoch12Epoch2(2000.0d, epoch2);
                result = m4.multiply(m3).multiply(m2).multiply(m1);               
                break;
            case FK4:            
            case FK4_NO_E:
                m1 = convertIAU2006MatrixEpoch12Epoch2(epoch1, 2000.0d);
                m2 = convertICRS2J2000Matrix().transpose();
                m3 = convertICRS2FK5Matrix();
                m4 = convertFK52FK4Matrix(epobs);
                m5 = convertBesselianMatrixEpoch12Epoch2(1950d, epoch2);
                result = m5.multiply(m4).multiply(m3).multiply(m2).multiply(m1);                
                break;
            case J2000:
                result = convertIAU2006MatrixEpoch12Epoch2(epoch1, epoch2);
                break;
            default:
                throw new JWcsError("J2000 to "+s2.getName()+" conversion is not supported");
        }
        return result;        
    }    

    /**
     * Converts FK4 to another coordinate reference frame.  
     * 
     * <p>The matrix to compute FK4 to <code>s2</code> is computed in this way:
     * <ul>
     * <li>for ICRS :  m3.multiply(m2).multiply(m1) with m1 = convertBesselianMatrixEpoch12Epoch2(epoch1, 1950.0),
     * m2 = convertFK42FK5Matrix(epobs) and m3 = convertICRS2FK5Matrix().transpose()</li>
     * <li>for FK5 : m3.multiply(m2).multiply(m1) with m1 = convertBesselianMatrixEpoch12Epoch2(epoch1, 1950.0),
     * m2 = convertFK42FK5Matrix(epobs) and m3 = convertJulianMatrixEpoch12Epoch2(2000.0, epoch2)</li>
     * <li>for FK4 and FK4-NO-E : convertBesselianMatrixEpoch12Epoch2(epoch1, epoch2)</li>
     * <li>for J2000 : m5.multiply(m4).multiply(m3).multiply(m2).multiply(m1)
     * with m1 = convertBesselianMatrixEpoch12Epoch2(epoch1, 1950.0), 
     * m2 = convertFK52FK4Matrix(epobs).transpose(), m3 = convertICRS2FK5Matrix().transpose(), 
     * m4 = convertICRS2J2000Matrix() and m5 = convertIAU2006MatrixEpoch12Epoch2(2000.0d, epoch2)</li>
     * </ul>          
     *
     * @param s2 Output coordinate reference frame
     * @param epoch1 Epoch belonging to coordinate reference frame s1 : either
     * Besselian or Julian.
     * @param epoch2 Epoch belonging to coordinate reference frame s2 : either 
     * Besselian or Julian.
     * @param epobs Epoch of observation. Only valid for conversions between FK4
     * and FK5.
     * @return Rotation matrix to transform a position in FK4 with <code>epoch1</code>
     * to an equatorial crs with equator and equinox at <code>epoch2</code> in 
     *  coordinate reference frame <code>s2</code>.
     * @throws JWcsError FK4 to <code>s2</code> conversion is not supported
     * @see #convertBesselianMatrixEpoch12Epoch2
     * @see #convertFK42FK5Matrix
     * @see #convertICRS2FK5Matrix
     * @see #convertJulianMatrixEpoch12Epoch2
     * @see #convertFK52FK4Matrix
     * @see #convertICRS2J2000Matrix
     * @see #convertIAU2006MatrixEpoch12Epoch2
     */    
    private static RealMatrix convertFK4To(final CoordinateReferenceFrame.ReferenceFrame s2, final double epoch1, final double epoch2, final double epobs) {
        final RealMatrix result;       
        final RealMatrix m1;
        final RealMatrix m2;
        final RealMatrix m3;        
        final RealMatrix m4;         
        final RealMatrix m5;         
        switch(s2) {
            case ICRS:
                m1 = convertBesselianMatrixEpoch12Epoch2(epoch1, 1950.0d);
                m2 = convertFK42FK5Matrix(epobs);
                m3 = convertICRS2FK5Matrix().transpose();
                result = m3.multiply(m2).multiply(m1);                             
                break;
            case FK5:
                m1 = convertBesselianMatrixEpoch12Epoch2(epoch1, 1950.0d);
                m2 = convertFK42FK5Matrix(epobs);
                m3 = convertJulianMatrixEpoch12Epoch2(2000.0d, epoch2);
                result = m3.multiply(m2).multiply(m1);            
                break;
            case FK4:            
            case FK4_NO_E:
                result = convertBesselianMatrixEpoch12Epoch2(epoch1, epoch2);               
                break;
            case J2000:
                m1 = convertBesselianMatrixEpoch12Epoch2(epoch1, 1950.0d);
                m2 = convertFK52FK4Matrix(epobs).transpose();
                m3 = convertICRS2FK5Matrix().transpose();
                m4 = convertICRS2J2000Matrix();
                m5 = convertIAU2006MatrixEpoch12Epoch2(2000.0d, epoch2);
                result = m5.multiply(m4).multiply(m3).multiply(m2).multiply(m1);              
                break;
            default:
                throw new JWcsError("FK4 to "+s2.getName()+" conversion is not supported");
        }
        return result;        
    }
    
    /**
     * Converts FK5 to another coordinate reference frame.     
     *
     * <p>The matrix to compute FK5 to <code>s2</code> is computed in this way:
     * <ul>
     * <li>for ICRS :  m2.multiply(m1) with m1 = convertJulianMatrixEpoch12Epoch2(epoch1, 2000.0),
     * m2 = convertICRS2FK5Matrix().transpose()</li>
     * <li>for FK5 : convertJulianMatrixEpoch12Epoch2(epoch1, epoch2)</li>
     * <li>for FK4 and FK4-NO-E : m3.multiply(m2).multiply(m1) with 
     * m1 = convertJulianMatrixEpoch12Epoch2(epoch1, 2000), m2 = convertFK52FK4Matrix(epobs),
     * m3 = convertBesselianMatrixEpoch12Epoch2(1950.0, epoch2)</li>
     * <li>for J2000 : m4.multiply(m3).multiply(m2).multiply(m1) with
     * m1 = convertJulianMatrixEpoch12Epoch2(epoch1, 2000.0),
     * m2 = convertICRS2FK5Matrix().transpose(), m3 = convertICRS2J2000Matrix(), 
     * m4 = convertIAU2006MatrixEpoch12Epoch2(2000.0, epoch2)</li>
     * </ul>   
     * 
     * @param s2 Output coordinate reference frame
     * @param epoch1 Epoch belonging to coordinate reference frame s1 : either 
     * Besselian or Julian.
     * @param epoch2 Epoch belonging to coordinate reference frame s2 : either 
     * Besselian or Julian.
     * @param epobs Epoch of observation. Only valid for conversions between FK4
     * and FK5.
     * @return Rotation matrix to transform a position in FK5 with <code>epoch1</code>
     * to an equatorial crs with equator and equinox at <code>epoch2</code> in 
     * coordinate reference frame <code>s2</code>.
     * @throws JWcsError FK5 to <code>s2</code> conversion is not supported
     * @see #convertJulianMatrixEpoch12Epoch2
     * @see #convertICRS2FK5Matrix
     * @see #convertFK52FK4Matrix
     * @see #convertBesselianMatrixEpoch12Epoch2
     * @see #convertICRS2J2000Matrix
     * @see #convertIAU2006MatrixEpoch12Epoch2
     */    
    private static RealMatrix convertFK5To(final CoordinateReferenceFrame.ReferenceFrame s2, final double epoch1, final double epoch2, final double epobs) {
        final RealMatrix result;       
        final RealMatrix m1;
        final RealMatrix m2;
        final RealMatrix m3;        
        final RealMatrix m4;         
        switch(s2) {
            case ICRS:
                m1 = convertJulianMatrixEpoch12Epoch2(epoch1, 2000.0d);
                m2 = convertICRS2FK5Matrix().transpose();
                result = m2.multiply(m1);                                
                break;
            case FK5:
                result = convertJulianMatrixEpoch12Epoch2(epoch1, epoch2);              
                break;
            case FK4:            
            case FK4_NO_E:
                m1 = convertJulianMatrixEpoch12Epoch2(epoch1, 2000d);
                m2 = convertFK52FK4Matrix(epobs);
                m3 = convertBesselianMatrixEpoch12Epoch2(1950.0d, epoch2);
                result = m3.multiply(m2).multiply(m1);                 
                break;
            case J2000:
                m1 = convertJulianMatrixEpoch12Epoch2(epoch1, 2000.0d);
                m2 = convertICRS2FK5Matrix().transpose();
                m3 = convertICRS2J2000Matrix();
                m4 = convertIAU2006MatrixEpoch12Epoch2(2000.0d, epoch2);
                result = m4.multiply(m3).multiply(m2).multiply(m1);              
                break;
            default:
                throw new JWcsError("FK5 to "+s2.getName()+" conversion is not supported");
        }
        return result;        
    }
    
    /**
     * Converts ICRS to another coordinate reference frame.    
     * 
     * <p>The matrix to compute ICRS to <code>s2</code> is computed in this way:
     * <ul>
     * <li>for ICRS :  createRealIdentityMatrix(3)</li>
     * <li>for FK5 : m2.multiply(m1) with m1 = convertICRS2FK5Matrix() and 
     * m2 = convertJulianMatrixEpoch12Epoch2(2000.0, epoch2)</li>
     * <li>for FK4 and FK4-NO-E : m3.multiply(m2).multiply(m1) with 
     * m1 = convertICRS2FK5Matrix(), m2 = convertFK52FK4Matrix(epobs) and 
     * m3 = convertBesselianMatrixEpoch12Epoch2(1950.0, epoch2)</li>
     * <li>for J2000 : m2.multiply(m1) with m1 = convertICRS2J2000Matrix() 
     * and m2 = convertIAU2006MatrixEpoch12Epoch2(2000.0d, epoch2)</li>
     * </ul>      
     *
     * @param s2 Output coordinate reference frame
     * @param epoch2 Epoch belonging to coordinate reference frame s2 : either 
     * Besselian or Julian.
     * @param epobs Epoch of observation. Only valid for conversions between FK4
     * and FK5.
     * @return Rotation matrix to transform a position in FK4 with <code>epoch1</code>
     * to an equatorial crs with equator and equinox at <code>epoch2</code> in 
     *  coordinate reference frame <code>s2</code>.
     * @throws JWcsError FK4 to <code>s2</code> conversion is not supported
     * @see NumericalUtility#createRealIdentityMatrix
     * @see #convertICRS2FK5Matrix
     * @see #convertJulianMatrixEpoch12Epoch2
     * @see #convertFK52FK4Matrix
     * @see #convertBesselianMatrixEpoch12Epoch2
     * @see #convertICRS2J2000Matrix
     * @see #convertIAU2006MatrixEpoch12Epoch2
     */    
    private static RealMatrix convertICRSTo(final CoordinateReferenceFrame.ReferenceFrame s2, final double epoch2, final double epobs) {
        final RealMatrix result;       
        final RealMatrix m1;
        final RealMatrix m2;
        final RealMatrix m3;        
        switch(s2) {
            case ICRS:
                result = createRealIdentityMatrix(3);
                break;
            case FK5:
                m1 = convertICRS2FK5Matrix();
                m2 = convertJulianMatrixEpoch12Epoch2(2000.0d, epoch2);
                result = m2.multiply(m1);                
                break;
            case FK4:            
            case FK4_NO_E:
                m1 = convertICRS2FK5Matrix();
                m2 = convertFK52FK4Matrix(epobs);
                m3 = convertBesselianMatrixEpoch12Epoch2(1950.0d, epoch2);
                result = m3.multiply(m2).multiply(m1);                
                break;
            case J2000:
                m1 = convertICRS2J2000Matrix();
                m2 = convertIAU2006MatrixEpoch12Epoch2(2000.0d, epoch2);
                result = m2.multiply(m1);                
                break;
            default:
                throw new JWcsError("ICRS to "+s2.getName()+" conversion is not supported");
        }
        return result;
    }
    
    /**
     * Remove the elliptic component of annual aberration when this is included
     * in a catalogue fk4 position.
     *
     * @param xyz vector xyz
     * @param eterm E-terms vector (as returned by getEterms()). If input is
     * omitted (== <code>null</code>), the e-terms for 1950 will be substituted.
     * @return Mean place
     * @see FK4#getEterms
     */
    private static RealMatrix removeEterms(final RealMatrix xyz, final RealMatrix eterm) {
        RealMatrix etermToRemove = eterm;
        if (eterm == null) {
            etermToRemove = FK4.getEterms(1950);
        }
        return xyz.subtract(etermToRemove.transpose());
    }

    /**
     * Add the elliptic component of annual aberration when the rsult must be a
     * catalogue fk4 position.    
     *
     * @param xyz Cartesian position(s) converted from long/lat
     * @param eterm E-terms vector (as returned by getEterms()). If input is
     * omitted (i.e. == <code>null</code>), the e-terms for 1950 will be substituted.
     * @return Apparent place
     * @see FK4#getEterms
     */
    private static RealMatrix addEterms(final RealMatrix xyz, final RealMatrix eterm) {
        RealMatrix etermToAdd = eterm;
        if (eterm == null) {
            etermToAdd = FK4.getEterms(1950);
        }
        return xyz.add(etermToAdd.transpose());
    }

    /**
     * Create matrix to convert equatorial fk4 coordinates (without e-terms) to
     * IAU 1958 lII,bII system of galactic coordinates.     
     *
     * <p>Original definitions from 1.:
     * <ul>
     * <li>The new north galactic pole lies in the direction 
     * alpha = 12h49m (192.25 deg), delta=27.4 deg (equinox 1950.0). </li>
     * <li>The new zero of longitude is the great semicircle
     * originating at the new north galactic pole at the position angle theta =
     * 123 deg with respect to the equatorial pole for 1950.0.</li>
     * <li>Longitude increases from 0 to 360 deg. The sense is such that,
     * on the galactic equator increasing galactic longitude corresponds to 
     * increasing Right Ascension. Latitude increases from -90 deg through 0 deg
     * to 90 deg at the new galactic pole.</li>
     * </ul>
     *
     * <p>Given the RA and Dec of the galactic pole, and using the Euler angles
     * scheme: M = rotZ(a3).rotY(a2).rotZ(a1)
     *
     * <p>We first rotate the spin vector of the XY plane about an angle a1 =
     * ra_pole and then rotate the spin vector in the XZ plane (i.e. around the
     * Y axis) with an angle a2=90-dec_pole to point it in the right
     * declination.
     *
     * <p>Now think of a circle with the galactic pole as its center. The radius is
     * equal to the distance between this center and the equatorial pole. The
     * zero point now is on the circle and opposite to this pole.
     *
     * <p>We need to rotate along this circle (i.e. a rotation around the new
     * z-axis) in a way that the angle between the zero point and the equatorial
     * pole is equal to 123 deg.
     * 
     * <p>So first we need to compensate for the 180 deg of the current zero
     * longitude, opposite to the pole. Then we need to rotate about an angle
     * 123 deg but in a way that increasing galactic longitude corresponds to
     * increasing Right Ascension which is opposite to the standard rotation of
     * this circle (note that we rotated the original X axis about 192.25 deg).
     * The last rotation angle therefore is a3=+180-123: M =
     * rotZ(180-123.0)*rotY(90-27.4)*rotZ(192.25)
     * 
     * <p>Reference:<br> 
     * ----------<br> 
     * 1. Blaauw, A., Gum C.S., Pawsey, J.L., Westerhout,
     * G.: 1958,<br> 
     * 2. Monthly Notices Roy. Astron. Soc. 121, 123, 3. Blaauw, A.,
     * 2007. Private communications.<br>     
     *
     * <p>The composed rotation matrix is the same as in Slalib's 'ge50.f' and the
     * matrix in eq. (32) of Murray (1989).
     *
     * @return 3x3 RealMatrix M as in XYZgal = M * XYZb1950
     */
    protected static final RealMatrix convertMatrixEqB19502Gal() {
        return rotZ(180d - 123.0d).multiply(rotY(90d - 27.4d)).multiply(rotZ(192.25d));
    }

    /**
     * Transforms galactic to supergalactic coordinates.     
     *
     * <p>The Supergalactic equator is conceptually defined by the
     * plane of the local (Virgo-Hydra-Centaurus) supercluster, and the origin
     * of supergalactic longitude is at the intersection of the supergalactic
     * and galactic planes. (de Vaucouleurs)
     *
     * <p>North SG pole at l=47.37 deg, b=6.32 deg. Node at l=137.37, sgl=0
     * (inclination 83.68 deg).
     *
     * <p>Older references give for he position of the SG node 137.29 which differs
     * from 137.37 deg in the official definition.
     *
     * <p>For the rotation matrix we chose the scheme <code>Rz.Ry.Rz</code> 
     * Then first we rotate about 47.37 degrees along the z-axis allowed by a 
     * rotation about 90-6.32 degrees is needed to set the pole to the right 
     * declination. The new plane intersects the old one at two positions. 
     * One of them is l=137.37, b=0 (in galactic coordinates). If we want this 
     * to be sgl=0 we have to rotate this plane along the new z-axis about an 
     * angle of 90 degrees. So the composed rotation matrix is: 
     * M = Rotz(90)*Roty(90-6.32)*Rotz(47.37) 
     * 
     * <p>Reference: <br> 
     * ----------<br> 
     * Lahav, O., The supergalactic plane revisited with
     * the Optical Redshift Survey Mon. Not. R. Astron. Soc. 312, 166-176 (2000)     
     *
     * @return RealMatrix M as in <code>XYZsgal = M * XYZgal</code>
     */
    protected static final RealMatrix convertMatrixGal2Sgal() {
        return rotZ(90.0d).multiply(rotY(90d - 6.32d)).multiply(rotZ(47.37d));
    }

    /**
     * What is the obliquity of the ecliptic at this Julian date? (IAU model
     * 2000).
     *
     * <p>The epoch is entered in Julian date and the time is
     * calculated w.r.t. J2000. The obliquity is the angle between the mean
     * equator and ecliptic, or, between the ecliptic pole and mean celestial
     * pole of date.<br>
     * <br>
     * Reference:<br> 
     * ----------<br>
     * Fukushima, T. 2003, AJ, 126,1 Kaplan, H., 2005, The
     * IAU Resolutions on Astronomical Reference Systems, TimeUtility Scales, and
     * Earth Rotation Models, United States Naval Observatory circular no. 179,
     * http://aa.usno.navy.mil/publications/docs/Circular_179.pdf (page 44)    
     *
     * @param jd Julian date
     * @return Mean obliquity in degrees
     */
    private static double obliquity2000(final double jd) {
        // T = (Date - 1 jan, 2000, 12h noon)
        final double T = (jd - 2451545.0d) / 36525.0d;
        final double eps = (84381.406d
                + (-46.836769d
                + (-0.0001831d
                + (0.00200340d
                + (-0.000000576d
                + (-0.0000000434d) * T) * T) * T) * T) * T) / 3600.0d;
        return eps;
    }

    /**
     * What is the obliquity of the ecliptic at this Julian date? (IAU 1980
     * model).
     *
     * <p>The epoch is entered in Julian date and the time is
     * calculated w.r.t. J2000. The obliquity is the angle between the mean
     * equator and ecliptic, or, between the ecliptic pole and mean celestial
     * pole of date<br>
     * <br>
     * Reference:<br>
     * ---------- <br>
     * Explanatory Supplement to the Astronomical Almanac,
     * P. Kenneth Seidelmann (ed), University Science Books (1992), Expression
     * 3.222-1 (p114).     
     *
     * @param jd Julian date
     * @return Mean obliquity in degrees
     */
    private static double obliquity1980(final double jd) {
        // T = (Date - 1 jan, 2000, 12h noon)
        final double T = (jd - 2451545.0d) / 36525.0d;
        return (84381.448d + (-46.8150d + (-0.00059d + 0.001813d * T) * T) * T) / 3600.0d;
    }

    /**
     * Calculates a rotation matrix to convert equatorial coordinates to
     * ecliptical coordinates.     
     * 
     * <p>1. The origin for ecliptic longitude is the vernal equinox.
     * Therefore the coordinates of a fixed object is subject to shifts due to
     * precession. The rotation matrix uses the obliquity to do the conversion
     * to the wanted ecliptic coordinates. So we always need to enter an epoch.
     * Usually this is J2000, but it can also be the epoch of date. The
     * additional reference system indicates whether we need a Besselian or a
     * Julian epoch.
     * <br>
     * 2. In the FITS paper of Calabretta and Greisen (2002), one observes the
     * following relations to FITS: 
     * <ul>
     * <li>Keyword RADESYSa sets the catalog system
     * FK4, FK4-NO-E or FK5 This applies to equatorial and ecliptical
     * coordinates with the exception of FK4-NO-E.</li>
     * <li>FK4 coordinates are not strictly spherical since they include a 
     * contribution from the elliptic terms of aberration, the so-called e-terms
     * which amount to max. 343 milliarcsec. FITS paper: <i>'Strictly speaking, 
     * therefore, a map obtained from, say, a radio synthesis telescope, should 
     * be regarded as FK4-NO-E unless it has been appropriately re-sampled or a 
     * distortion correction provided. In common usage, however, CRVALia for 
     * such maps is usually given in FK4 coordinates. In doing so, the e-terms 
     * are effectively corrected to first order only.'</i>. 
     * (See also ES, eq. 3.531-1 page 170).</li>
     * <li>Keyword EQUINOX sets the epoch of the mean equator and equinox.</li> 
     * <li>Keyword EPOCH is often used in older FITS files. It is a deprecated 
     * keyword and should be replaced by EQUINOX. It does not require keyword 
     * RADESYS. From its value we derive whether the reference system is FK4 or FK5 (the
     * marker value is 1984.0)</li>
     * <li>Ecliptic coordinates require the epoch of the
     * equator and equinox of date. This will be taken as the time of
     * observation rather than EQUINOX. FITS paper: <i>'The time of observation
     * may also be required for other astrometric purposes in addition to the
     * usual astrophysical uses, for example, to specify when the mean place was
     * correct in accounting for proper motion, including "fictitious" proper
     * motions in the conversion between the FK4 and FK5 systems. The old
     * *DATE-OBS* keyword may be used for this purpose. However, to provide a
     * more convenient specification we here introduce the new keyword
     * MJD-OBS'.</i> So MJD-OBS is the modified Julian Date (JD - 2400000.5) of 
     * the start of the observation.</li>
     * </ul>
     *<br>
     * 3. Equatorial to ecliptic transformations use the time dependent
     * obliquity of the equator (also known as the obliquity of the ecliptic).
     * Again, start with: <code>M = rotZ(0).rotX(eps).rotZ(0) = E.rotX(eps).E =
     * rotX(eps)</code> In fact this is only a rotation around the X axis
     * 
     * <p>Reference:<br>
     * ----------<br>
     * Representations of celestial coordinates in FITS, 
     * Calabretta. M.R., and Greisen, E.w., (2002) Astronomy and Astrophysics,
     * 395, 1077-1122. http://www.atnf.csiro.au/people/mcalabre/WCS/ccs.pdf     
     *
     * @param epoch Epoch of the equator and equinox of date
     * @param refSystem equatorial refrence frame to determine if one entered epoch in B
     * or J coordinates
     * @return 3x3 RealMatrix M as in <code>XYZecl = M * XYZeq</code>
     */
    protected static final RealMatrix convertMatrixEq2Ecl(final double epoch, final CoordinateReferenceFrame.ReferenceFrame refSystem) {
        final double jd;
        if (CoordinateReferenceFrame.ReferenceFrame.FK4.equals(refSystem)) {
            jd = convertEpochBessel2JD(epoch);
        } else {
            jd = convertEpochJulian2JD(epoch);
        }
        final double eps;
        if (CoordinateReferenceFrame.ReferenceFrame.ICRS.equals(refSystem) || CoordinateReferenceFrame.ReferenceFrame.J2000.equals(refSystem)) {
            eps = obliquity2000(jd);
        } else {
            eps = obliquity1980(jd);
        }
        return rotX(eps);
    }  
    
    /**
     * Precession from one epoch to another in the fk5 coordinate reference frame.
     *
     * <p>The precession matrix is: <code>M = rotZ(-z).rotY(\u03B8).rotZ(-\u03B6)</code><br>
     * <br> 
     * Reference:<br> 
     * ----------<br> 
     * Seidelman, P.K., 1992. Explanatory Supplement to
     * the Astronomical Almanac. University Science Books, Mill Valley. 3.214 p
     * 106
     *
     * @param jEpoch1 Julian start epoch
     * @param jEpoch2 Julian epoch to process to
     * @return 3x3 rotation matrix M as in <code>XYZepoch2 = M * XYZepoch1</code>
     */
    private static RealMatrix convertJulianMatrixEpoch12Epoch2(final double jEpoch1, final double jEpoch2) {
        final double jd1 = convertEpochJulian2JD(jEpoch1);
        final double jd2 = convertEpochJulian2JD(jEpoch2);
        final double[] precessionAngles = lieskeprecangles(jd1, jd2);
        return precessionMatrix(precessionAngles[0], precessionAngles[1], precessionAngles[2]);
    }

    /**
     * Calculates IAU 1976 precession angles for a precession of epoch
     * corresponding to Julian date jd1 to epoch corresponds to Julian date jd2.
     *
     * <p>References:<br>
     * Lieske,J.H., 1979. Astron.Astrophys.,73,282.<br>
     * equations (6) and (7), p283.<br>
     * Kaplan,G.H., 1981. USNO circular no. 163, pA2.
     *
     * @param jd1 Julian date for start epoch
     * @param jd2 Julian date for end epoch
     * @return Angles \u03B6 (zeta), z, \u03B8 (theta) in decimal degrees
     */
    private static double[] lieskeprecangles(final double jd1, final double jd2) {
        // T = (Current epoch - 1 jan, 2000, 12h noon)
        final double T0 = (jd1 - 2451545.0d) / 36525.0d;
        final double T = (jd2 - jd1) / 36525.0d;

        final double w = 2306.2181d+(1.39656d-0.000139d*T0)*T0;
        final double zeta = (w+((0.30188d-0.000344d*T0)+0.017998d*T)*T)*T;
        final double z = (w+(1.09468d+0.000066d*T0+0.018203d*T)*T)*T;
        final double theta = (2004.3109d+(-0.85330d-0.000217d*T0)*T0+
                ((-0.42665d-0.000217d*T0)-0.041833d*T)*T)*T;        
        //Return values in degrees
        return new double[]{zeta / 3600.0d, z / 3600.0d, theta / 3600.0d};
    }

    /**
     * Given three precession angles, creates the corresponding rotation matrix.
     *
     * <p>Return the precession matrix for the three precession angles zeta, z and
     * theta. Rotation matrix: <code>M = rotZ(-z).rotY(\u03B8).rotZ(-\u03B6)</code> (ES 3.21-7, p
     * 103). Also allowed is the expression: <code>rotZ(-90-z)*rotX(th)*rotZ(90-zeta)</code>
     *
     * @param zeta \u03B6  in decimal degrees
     * @param z z in decimal degree
     * @param theta \u03B8 in decimal degrees
     * @return Rotation matrix M as in <code>XYZepoch1 = M * XYZepoch2</code>
     */
    private static RealMatrix precessionMatrix(final double zeta, final double z, final double theta) {
        return rotZ(-z).multiply(rotY(theta)).multiply(rotZ(-zeta));
    }

    /**
     * Precession from one epoch to another in the FK4 coordinate reference frame.
     *
     * <p>The precession matrix is: <code>M = rotZ(-z).rotY(\u03B8).rotZ(-\u03B6)</code><br>
     * <br>
     * Reference:<br> 
     * ----------<br> 
     * Seidelman, P.K., 1992. Explanatory Supplement to
     * the Astronomical Almanac. University Science Books, Mill Valley. 3.214 p
     * 106
     *
     * @param bEpoch1 Besselian start epoch
     * @param bEpoch2 Besselian epoch to precess to.
     * @return 3x3 rotation matrix M as in <code>XYZepoch2 = M * XYZepoch1</code>
     * @see #newcombPrecAngles
     * @see #precessionMatrix
     */
    private static RealMatrix convertBesselianMatrixEpoch12Epoch2(final double bEpoch1, final double bEpoch2) {
        final double[] precessionAngles = newcombPrecAngles(bEpoch1, bEpoch2);
        return precessionMatrix(precessionAngles[0], precessionAngles[1], precessionAngles[2]);
    }

    /**
     * Calculates precession angles for a precession in FK4, using Newcomb's
     * method (Woolard and Clemence angles).
     *
     * <p>Newcomb's precession angles for old catalogs (FK4), see ES
     * 3.214 p.106. Input are <b>Besselian epochs</b>! Adopted accumulated
     * precession angles from equator and equinox at B1950 to 1984 January 1d 0h
     * according to ES (table 3.214.1, p 107) are: <code>zeta=783.7092</code>, 
     * <code>z=783.8009</code> and <code>theta=681.3883</code>. The Woolard and 
     * Clemence angles (derived in this routine) are: <code>zeta=783.70925</code>,
     * <code>z=783.80093</code> and <code>theta=681.38830</code>
     *
     * @param epoch1 Besselian start epoch
     * @param epoch2 Besselian end epoch
     * @return Angles \u03B6 (zeta), z, \u03B8 (theta) degrees
     */
    private static double[] newcombPrecAngles(final double epoch1, final double epoch2) {
        final double t1 = (epoch1 - 1850.0d) / 1000.0d;    //1000 tropical years
        final double t2 = (epoch2 - 1850.0d) / 1000.0d;
        final double tau = t2 - t1;

        double d0 = 23035.545d;
        double d1 = 139.720d;
        double d2 = 0.060d;
        double d3 = 30.240d;
        double d4 = -0.27d;
        double d5 = 17.995d;
        double a0 = d0 + t1 * (d1 + d2 * t1);
        double a1 = d3 + d4 * t1;
        double a2 = d5;
        final double zeta_a = tau * (a0 + tau * (a1 + tau * a2));

        d0 = 23035.545d;
        d1 = 139.720d;
        d2 = 0.060d;
        d3 = 109.480d;
        d4 = 0.39d;
        d5 = 18.325d;
        a0 = d0 + t1 * (d1 + d2 * t1);
        a1 = d3 + d4 * t1;
        a2 = d5;
        final double z_a = tau * (a0 + tau * (a1 + tau * a2));

        d0 = 20051.12d;
        d1 = -85.29d;
        d2 = -0.37d;
        d3 = -42.65d;
        d4 = -0.37d;
        d5 = -41.80d;
        a0 = d0 + t1 * (d1 + d2 * t1);
        a1 = d3 + d4 * t1;
        a2 = d5;
        final double theta_a = tau * (a0 + tau * (a1 + tau * a2));
        // Return values in degrees
        return new double[]{zeta_a / 3600.0d, z_a / 3600.0d, theta_a / 3600.0d};
    }

    /**
     * Creates a matrix to precess from B1950 in FK4 to J2000 in FK5 following
     * to Murray's (1989) procedure.
     *
     * <p>Murray precesses from B1950 to J2000 using a precession
     * matrix by Lieske. Then applies the equinox correction and ends up with a
     * transformation matrix <b>X(0)</b> as given in this function. In Murray's
     * article it is proven that using the procedure as described in the
     * article, <code>r_fk5 = X(0).r_fk4</code> for extra galactic sources where we
     * assumed that the proper motion in FK5 is zero. This procedure is
     * independent of the epoch of observation. Note that the matrix is not a
     * rotation matrix. FK4 is not an inertial coordinate frame (because of the
     * error in precession and the motion of the equinox. This has consequences
     * for the proper motions. e.g. a source with zero proper motion in FK5 has
     * a fictitious proper motion in FK4. This affects the actual positions in a
     * way that the correction is bigger if the epoch of observation is further
     * away from 1950.0 The focus of this library is on data of which we do not
     * have information about the proper motions. So for positions of which we
     * allow non zero proper motion in FK5 one needs to supply the epoch of
     * observation.<br>
     * <br>
     * Reference:<br> 
     * ----------<br>
     * <ul>
     * <li>Murray, C.A. The Transformation of coordinates
     * between the systems B1950.0 and J2000.0, and the principal galactic axis
     * referred to J2000.0, Astronomy and Astrophysics (ISSN 0004-6361), vol.
     * 218, no. 1-2, July 1989, p. 325-329.</li> 
     * <li>Poppe P.C.R.,, Martin, V.A.F.,
     * Sobre as Bases de Referencia Celeste SitientibusSerie Ciencias Fisicas</li>
     * </ul>     
     *
     * @param t Besselian epoch as epoch of observation
     * @return 3x3 matrix M as in <i>XYZfk5 = M * XYZfk4</i>
     * @see TimeUtility#convertEpochBessel2JD
     */
    private static RealMatrix convertFK42FK5Matrix(final Double t) {
        final RealMatrix mat = convertFK42FK5Matrix();
        if (!Double.isNaN(t)) {
            final double jd = convertEpochBessel2JD(t);
            final double T = (jd - 2433282.423d) / 36525.0d; //t-1950 in Julian centuries = F^-1.t1 from Murray (1989)
            final double r00 = mat.getEntry(0, 0) - 0.0026455262d * T / 1000000.0d;
            final double r01 = mat.getEntry(0, 1) - 1.1539918689d * T / 1000000.0d;
            final double r02 = mat.getEntry(0, 2) + 2.1111346190d * T / 1000000.0d;
            final double r10 = mat.getEntry(1, 0) + 1.1540628161d * T / 1000000.0d;
            final double r11 = mat.getEntry(1, 1) - 0.0129042997d * T / 1000000.0d;
            final double r12 = mat.getEntry(1, 2) + 0.0236021478d * T / 1000000.0d;
            final double r20 = mat.getEntry(2, 0) - 2.1112979048d * T / 1000000.0d;
            final double r21 = mat.getEntry(2, 1) - 0.0056024448d * T / 1000000.0d;
            final double r22 = mat.getEntry(2, 2) + 0.0102587734d * T / 1000000.0d;
            mat.setEntry(0, 0, r00);
            mat.setEntry(0, 1, r01);
            mat.setEntry(0, 2, r02);
            mat.setEntry(1, 0, r10);
            mat.setEntry(1, 1, r11);
            mat.setEntry(1, 2, r12);
            mat.setEntry(2, 0, r20);
            mat.setEntry(2, 1, r21);
            mat.setEntry(2, 2, r22);
        }
        return mat;
    }

    /**
     * Create a matrix to precess from B1950 in FK4 to J2000 in FK5 following to
     * Murray’s (1989) procedure without epoch of observation.
     *
     * @return 3x3 matrix M as in <code>XYZfk5 = M * XYZfk4</code>
     */
    private static RealMatrix convertFK42FK5Matrix() {
        final double[][] array = {
            {0.9999256794956877d, -0.0111814832204662d, -0.0048590038153592d},
            {0.0111814832391717d, 0.9999374848933135d, -0.0000271625947142d},
            {0.0048590037723143d, -0.0000271702937440d, 0.9999881946023742d}
        };
        return createRealMatrix(array);
    }

    /**
     * Creates a matrix to convert a position in fk5 to fk4 using the inverse 
     * matrix convertFK42FK5Matrix.
     *
     * @param t Epoch of observation for those situations where we allow no-zero
     * proper motion in fk4
     * @return Rotation matrix M as in <code>XYZfk5 = M * XYZfk4</code>
     */
    private static RealMatrix convertFK52FK4Matrix(final Double t) {
        return inverse(convertFK42FK5Matrix(t));
    }

    /**
     * Creates a rotation matrix to convert a position from ICRS to fk5, J2000.
     *
     * <p>Notes:<br>
     * ------<br>
     * Return a matrix that converts a position vector in ICRS to
     * FK5, J2000. We do not use the first or second order approximations given
     * in the reference, but use the three rotation matrices from the same paper
     * to obtain the exact result <code>M = rotX(-eta0)*rotY(xi0)*rotZ(da0)</code>
     * <code>eta0 = -19.9 mas</code>, <code>xi0 = 9.1 mas</code> and 
     * <code>da0 = -22.9 mas</code><br>
     * <br>
     * Reference:<br> 
     * ----------<br>
     * Kaplan G.H., The IAU Resolutions on Astronomical
     * Reference systems, TimeUtility scales, and Earth Rotation Models, US Naval
     * Observatory, Circular No. 179     
     *
     * @return 3x3 rotation matrix M as in <code>XYZfk5 = M * XYZicrs</code>
     */
    private static RealMatrix convertICRS2FK5Matrix() {
        final double eta0 = -19.9d / (3600d * 1000d);  //Convert mas to degree
        final double xi0 = 9.1d / (3600d * 1000d);
        final double da0 = -22.9d / (3600d * 1000d);
        return rotX(-eta0).multiply(rotY(xi0)).multiply(rotZ(da0));
    }

    /**
     * Returns a rotation matrix for conversion of a position in the ICRS to the
     * dynamical reference system based on the dynamical mean equator and
     * equinox of J2000.0 (called the dynamical J2000 system).
     *
     * <p>Return a matrix that converts a position vector in ICRS to
     * Dyn. J2000. We do not use the first or second order approximations given
     * in the reference, but use the three rotation matrices to obtain the
     * exact.<br> 
     * <br>
     * Reference:<br>
     * ----------<br>
     * Hilton and Hohenkerk (2004), Astronomy and
     * Astrophysics 413, 765-770 Kaplan G.H., The IAU Resolutions on
     * Astronomical Reference systems, TimeUtility scales, and Earth Rotation
     * Models, US Naval Observatory, Circular No. 179
     * 
     * <p>Note the precession is applied only to equatorial
     * coordinates in the system of dynamical J2000 coordinates. When converting
     * from ICRS coordinates this means applying a frame bias. Therefore the
     * angles differ from the precession Fukushima-Williams angles (IAU 2006).
     * The precession matrix is: <code>M = rotZ(-z).rotY(\u03B8).rotZ(\u03B6)</code><br>
     * <br>
     * Reference:<br>
     * ----------<br>
     * Capitaine N. et al.: IAU 2000 precession A and A 412, 567-586 (2003)
     *
     * @return Rotation matrix to transform positions from ICRS to dyn J2000
     */
    private static RealMatrix convertICRS2J2000Matrix() {
        final double eta0 = -6.8192d / (3600d * 1000d); //Convert mas to degree
        final double xi0 = -16.617d / (3600d * 1000d);
        final double da0 = -14.6d / (3600d * 1000d);
        return rotX(-eta0).multiply(rotY(xi0)).multiply(rotZ(da0));
    }

    /**
     * Creates a rotation matrix for a precession based on IAU 2000/2006
     * expressions, see {@link #convertIAU2006PrecAngles(double) }.
     *
     * <p>Note that we apply this precession only to equatorial
     * coordinates in the system of dynamical J2000 coordinates. When converting
     * from ICRS coordinates this means applying a frame bias. Therefore the
     * angles differ from the precession Fukushima-Williams angles (IAU 2006)
     * The precession matrix is: <code>M = rotZ(-z).rotY(\u03B8).rotZ(-\u03B6)</code>
     * 
     * <p>Reference:<br>
     * ----------<br>
     * Capitaine N. et al.: IAU 2000 precession A and A 412, 567-586 (2003)
     *
     * @param epoch1 Julian start epoch
     * @param epoch2 Julian epoch to precess to
     * @return RealMatrix to transform equatorial coordinates from epoch1 to epoch2
     * as in <code>XYZepoch2 = M * XYZepoch1</code>
     */
    private static RealMatrix convertIAU2006MatrixEpoch12Epoch2(final double epoch1, final double epoch2) {
        final RealMatrix result;
        if (equal(epoch1, epoch2)) {
            result = createRealIdentityMatrix(3);
        } else if (equal(epoch1,2000.0,1e-3)) {
            final double[] precessionAngles = convertIAU2006PrecAngles(epoch2);
            result = precessionMatrix(precessionAngles[0], precessionAngles[1], precessionAngles[2]);
        } else { // If both epochs are not J2000.0
            double[] precessionAngles = convertIAU2006PrecAngles(epoch1);
            RealMatrix m1 = precessionMatrix(precessionAngles[0], precessionAngles[1], precessionAngles[2]);
            m1 = m1.transpose();
            precessionAngles = convertIAU2006PrecAngles(epoch2);
            final RealMatrix m2 = precessionMatrix(precessionAngles[0], precessionAngles[1], precessionAngles[2]);
            result = m1.multiply(m2);
        }
        return result;
    }

    /**
     * Calculates IAU 2000 precession angles for precession from input epoch to
     * J2000.
     *
     * <p>Input are Julian epochs! <code>T = (jd-2451545.0)/36525.0</code> 
     * Combined with <code>jd = Jepoch-2000.0)*365.25 + 2451545.0</code> gives:
     * <code>T = (epoch-2000.0)/100.0</code>
     * 
     * <p>Reference:<br>
     * ----------<br>
     * Capitaine N. et al., IAU 2000 precession A and A 412, 567-586 (2003)     
     *
     * @param epoch Julian epoch of observation
     * @return Angles \u03B6 (zeta), z, \u03B8 (theta) in degrees to setup a
     * rotation matrix to transform from J2000 to input epoch.     
     */
    private static double[] convertIAU2006PrecAngles(final double epoch) {
        // T = (Current epoch - 1 jan, 2000, 12h noon)
        final double T = (epoch - 2000.0d) / 100.0d;
        double d0 = 2.5976176d;
        double d1 = 2306.0809506d;
        double d2 = 0.3019015d;
        double d3 = 0.0179663d;
        double d4 = -0.0000327d;
        double d5 = -0.0000002d;
        final double zeta_a = T * (d1 + T * (d2 + T * (d3 + T * (d4 + T * d5)))) + d0;

        d0 = -2.5976176d;
        d1 = 2306.0803226d;
        d2 = 1.0947790d;
        d3 = 0.0182273d;
        d4 = 0.0000470d;
        d5 = -0.0000003d;
        final double z_a = T * (d1 + T * (d2 + T * (d3 + T * (d4 + T * d5)))) + d0;

        d0 = 0.0;
        d1 = 2004.1917476d;
        d2 = -0.4269353d;
        d3 = -0.0418251d;
        d4 = -0.0000601d;
        d5 = -0.0000001d;
        final double theta_a = T * (d1 + T * (d2 + T * (d3 + T * (d4 + T * d5)))) + d0;

        //Return values in degrees
        return new double[]{zeta_a / 3600.0d, z_a / 3600.0d, theta_a / 3600.0d};
    }    
    
    /**
     * Given two angles in longitude and latitude returns corresponding
     * Cartesian coordinates x,y,z.
     *
     * <p>NotesThe three coordinate axes x, y and z, the set of
     * right-handed Cartesian axes that correspond to the usual celestial
     * spherical coordinate system. The xy-plane is the equator, the z-axis
     * points toward the north celestial pole, and the x-axis points toward the
     * origin of right ascension.
     *
     * @param longitude longitude in decimal degrees
     * @param latitude latitude in decimal degrees
     * @return Corresponding values of x,y,z in same order as input
     */
    public static RealMatrix longlat2xyz(final double longitude, final double latitude) {
        return longlatRad2xyz(FastMath.toRadians(longitude), FastMath.toRadians(latitude));
    }

    /**
     * Given two angles in longitude and latitude returns corresponding
     * Cartesian coordinates x,y,z.
     *
     * <p>The three coordinate axes x, y and z, the set of
     * right-handed Cartesian axes that correspond to the usual celestial
     * spherical coordinate system. The xy-plane is the equator, the z-axis
     * points toward the north celestial pole, and the x-axis points toward the
     * origin of right ascension.
     *
     * @param longitudeRad longitude in radians
     * @param latitudeRad latitude in radians
     * @return Corresponding values of x,y,z in same order as input
     */
    public static RealMatrix longlatRad2xyz(final double longitudeRad, final double latitudeRad) {
        final double x = FastMath.cos(longitudeRad) * FastMath.cos(latitudeRad);
        final double y = FastMath.sin(longitudeRad) * FastMath.cos(latitudeRad);
        final double z = FastMath.sin(latitudeRad);
        final double[][] array = {
            {x},
            {y},
            {z}
        };
        return createRealMatrix(array);
    }

    /**
     * Given Cartesian x,y,z return corresponding longitude and latitude in
     * degrees.
     *
     * <p>Note that one can expect strange behavior for the values of
     * the longitudes very close to the pole. In fact, at the poles itself, the
     * longitudes are meaningless.
     *
     * @param xyz Vector with values for x,y,z
     * @return The same number of positions (longitude, latitude and in the same
     * order as the input.
     */
    public static double[] xyz2longlat(final RealMatrix xyz) {
        final double[] vec = xyz.getColumn(0);
        final double len = FastMath.sqrt(FastMath.pow(vec[0], 2)+FastMath.pow(vec[1], 2)+FastMath.pow(vec[2], 2));
        final double x = vec[0]/len;
        final double y = vec[1]/len;
        final double z = vec[2]/len;
        double longitude = FastMath.toDegrees(aatan2(y, x, 0));
        longitude = longitude < 0 ? longitude + 360.0d : longitude;
        final double latitude = FastMath.toDegrees(aasin(z));
        return new double[]{longitude, latitude};
    }    

    @Override
    public String toString() {
        final String result;
        final CoordinateReferenceFrame refSystem = getCoordinateReferenceFrame();
        if (refSystem == null) {
            result = getCoordinateReferenceSystem().getName();
        } else {
            result = getCoordinateReferenceSystem().getName()+" with "+refSystem.getReferenceFrame();
        }
        return result;
    }    
}
