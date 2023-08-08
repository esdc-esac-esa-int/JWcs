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
package io.github.malapert.jwcs;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.URI;
import java.net.URISyntaxException;
import java.text.ParseException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import io.github.malapert.jwcs.crs.AbstractCrs;
import io.github.malapert.jwcs.crs.Ecliptic;
import io.github.malapert.jwcs.crs.Equatorial;
import io.github.malapert.jwcs.crs.Galactic;
import io.github.malapert.jwcs.datum.FK4;
import io.github.malapert.jwcs.datum.FK5;
import io.github.malapert.jwcs.datum.ICRS;
import io.github.malapert.jwcs.datum.J2000;
import io.github.malapert.jwcs.position.SkyPosition;
import io.github.malapert.jwcs.proj.exception.JWcsException;
import io.github.malapert.jwcs.utility.DMS;
import io.github.malapert.jwcs.utility.HMS;
import io.github.malapert.jwcs.utility.HeaderFitsReader;
import io.github.malapert.jwcs.utility.NumericalUtility;
import io.github.malapert.jwcs.utility.TimeUtility;
import nom.tam.fits.Fits;
import nom.tam.fits.Header;
import nom.tam.fits.HeaderCard;
import nom.tam.util.Cursor;

/**
 * Tools to perform common operations without knowing the internals of JWcs. This class is a 
 * very basic summary of what can be done with JWcs, other specific transformations can use 
 * a similar code respect the one implemented here.<BR>
 * 
 * The tools provided are: conversion between different angle units (not specific to JWcs), 
 * conversion between different coordinate systems, precession between two epochs for coordinates 
 * referred to a given coordinate frame, conversions between B1950 and ICRS or FK5 frames, 
 * angle separation between two sky positions, conversions between Julian day numbers and ISO 
 * dates, conversion from decimal angles as degrees or hours into the tupple (deg, arcmin, arcsec 
 * or hour, min, sec), and a code to read the header of a FITS file to get ready to convert 
 * between sky and pixel coordinates.
 *
 * @author T. Alonso Albi (tomas.alonsoalbi@esa.int)
 */
public class JWcsTools {
    
    private static final String B1950 = "B1950";
    private static final int FIFTEEN = 15;
    private static final int SIXTY = 60;
    private static final int ONE = 1;

    private JWcsTools() {}

    /**
     * Basic coordinate systems to convert positions in the sky.
     */
    public enum COORD_SYS {
        EQUATORIAL,
        ECLIPTIC,
        GALACTIC
    }
    
    /**
     * Basic coordinate frames.
     */
    public enum FRAME {
        FK4,
        FK5,
        J2000
    }
    
    /**
     * Angle units to convert between units or arc and time.
     */
    public enum ANGLE {
        DEGREES (FIFTEEN),
        ARCMINUTES ((double) FIFTEEN * SIXTY),
        ARCSECONDS ((double) FIFTEEN * SIXTY * SIXTY),
        HOURS (ONE),
        MINUTES (SIXTY),
        SECONDS ((double) SIXTY * SIXTY),
        RADIANS (Math.toRadians(FIFTEEN))
        ;
                
        private double relativeConversionFactor;
        
        private ANGLE(double c) {
            relativeConversionFactor = c;
        }
        
        /**
         * Return the conversion factor for this specific unit into another unit.
         * @param to Desired output unit.
         * @return Conversion factor to that unit.
         */
        public double getConversionFactorTo(ANGLE to) {
            return to.relativeConversionFactor / relativeConversionFactor;            
        }
    }
    
    /**
     * Converts an angle from some unit into another.
     * @param angle The angle value.
     * @param from The unit for the previous value.
     * @param to The output unit desired.
     * @return The output angle for the output unit selected.
     */
    public static double convertAngle(double angle, ANGLE from, ANGLE to) {
        return angle * from.getConversionFactorTo(to);
    }
    
    /**
     * Returns the {@linkplain SkyPosition} object for a given direction and coordinate system.
     * @param lon Longitude or Ra in degrees.
     * @param lat Latitude or Dec in degrees.
     * @param c Coordinate system of the input position.
     * @return Instance of {@linkplain SkyPosition}.
     */
    public static SkyPosition getSkyPosition(double lon, double lat, COORD_SYS c) {
        if (lon < 0 || lon > 360)
            lon = Math.toDegrees(NumericalUtility.normalizeLongitude(Math.toRadians(lon)));
        SkyPosition input = null;
        switch (c) {
        case EQUATORIAL:
            input = new SkyPosition(lon, lat, new Equatorial(new FK5()));
            break;
        case ECLIPTIC:
            input = new SkyPosition(lon, lat, new Ecliptic());
            break;
        case GALACTIC:
            input = new SkyPosition(lon, lat, new Galactic());
            break;
        }
        return input;
    }
    
    /**
     * Transforms a given position in some coordinate system into another coordinate system.
     * @param input Input sky position in some coordinate system.
     * @param to Coordinate system for the output position.
     * @return Output position in degrees.
     */
    public static SkyPosition convert(SkyPosition input, COORD_SYS to) {
        switch (to) {
        case EQUATORIAL:
            return AbstractCrs.convertTo(new Equatorial(new FK5()), input);
        case ECLIPTIC:
            return AbstractCrs.convertTo(new Ecliptic(), input);
        case GALACTIC:
            return AbstractCrs.convertTo(new Galactic(), input);
        }
        
        return null; // Will never happen
    }
    
    /**
     * Computes the sky separation between two sky positions, possible defined respect different 
     * coordinate systems. See also {@linkplain NumericalUtility#distAngle(double[], double[])}.
     * @param pos1 Position 1.
     * @param pos2 Position 2.
     * @return Separation in degrees.
     */
    public static double skySeparation(SkyPosition pos1, SkyPosition pos2) {
        return SkyPosition.separation(pos1, pos2);
    }

    /**
     * Precess coordinates between two epochs, following the adequate method for a given coordinate frame.
     * @param ra Ra in degrees.
     * @param dec Dec in degrees.
     * @param inputEpoch Input epoch as a year with decimals.
     * @param outputEpoch Output epoch as a year with decimals.
     * @param frame The Frame in which the conversion will be performed. Note FK4 will use a very old
     * precession algorithm adequate only for very old catalogs and surveys, FK5 will use the classic 
     * IAU1976 formulae, and J2000 will use IAU2006 formulae, adequate for modern catalogs and surveys.
     * @return
     */
    public static double[] precessEquinox(double ra, double dec, double inputEpoch, double outputEpoch, FRAME frame) {
        AbstractCrs input = null;
        AbstractCrs output = null;
        switch (frame) {
        case FK4:
            input = new Equatorial(new FK4(B1950, "B"+inputEpoch));
            output = new Equatorial(new FK4(B1950, "B"+outputEpoch));
            break;
        case FK5:
            input = new Equatorial(new FK5("J"+inputEpoch));
            output = new Equatorial(new FK5("J"+outputEpoch));
            break;
        case J2000:
            J2000 in = new J2000();
            J2000 out = new J2000();
            in.setEquinox("J"+inputEpoch);
            out.setEquinox("J"+outputEpoch);
            input = new Equatorial(in);
            output = new Equatorial(out);
            break;
        }
        return AbstractCrs.convertTo(output, new SkyPosition(ra, dec, input)).getDoubleArray();
    }
    
    /**
     * Transform B1950 coordinates (equinox and epoch B1950) to ICRS.
     * @param ra Ra in degrees.
     * @param dec Dec in degrees.
     * @return Output ICRS position in degrees.
     */
    public static double[] convertB1950ToICRS(double ra, double dec) {
        return convertBxxxToICRS(ra, dec, 1950);
    }

    /**
     * Transform B1950 coordinates with a given observation epoch to ICRS. The observation 
     * epoch in usually B1950 also, but in some cases not. For instance, for IRAS data the 
     * observation epoch is B1983.5.
     * @param ra Ra in degrees.
     * @param dec Dec in degrees.
     * @param epochObs Observation epoch as a year.
     * @return Output ICRS position in degrees.
     */
    public static double[] convertBxxxToICRS(double ra, double dec, double epochObs) {
        if (ra < 0 || ra > 360)
            ra = Math.toDegrees(NumericalUtility.normalizeLongitude(Math.toRadians(ra)));
        SkyPosition pBxxxx = new SkyPosition(ra, dec, new Equatorial(new FK4(B1950, "B"+epochObs)));
        SkyPosition pICRS = AbstractCrs.convertTo(new Equatorial(new ICRS()), pBxxxx);

        return new double[] {pICRS.getLongitude(), pICRS.getLatitude()};
    }

    /**
     * Transform B1950 (equinox and epoch B1950) coordinates to FK5 J2000.
     * @param ra Ra in degrees.
     * @param dec Dec in degrees.
     * @return Output FK5 position in degrees.
     */
    public static double[] convertB1950ToFK5(double ra, double dec) {
        if (ra < 0 || ra > 360)
            ra = Math.toDegrees(NumericalUtility.normalizeLongitude(Math.toRadians(ra)));
        SkyPosition pB1950 = new SkyPosition(ra, dec, new Equatorial(new FK4()));
        SkyPosition pJ2000 = AbstractCrs.convertTo(new Equatorial(new FK5()), pB1950);

        return new double[] {pJ2000.getLongitude(), pJ2000.getLatitude()};
    }
        
    /**
     * Transform a Julian date to ISO date, for instance 2456915 to 2014-09-14T12:00:00.
     * @param jd Julian day number.
     * @return ISO date.
     */
    public static String julianDateToISO(double jd) {
        return TimeUtility.convertJulianDateToISO(jd);
    }
    
    /**
     * Transform an ISO date to a Julian day, for instance 2014-09-14T12:00:00 to 2456915.
     * @param iso The ISO date.
     * @return Julian day number.
     */
    public static double isoToJulianDate(String iso) throws ParseException {
        return TimeUtility.convertISOToJulianDate(iso);
    }
    
    /**
     * Reads a fits file and returns the {@linkplain JWcsMap} object to perform 
     * pixel <--> sky coordinate conversions.
     * @param url Path to the fits as a url, with file:// or http://.
     * @param hdu HDU number.
     * @return The instance of {@linkplain JWcsMap}.
     * @throws JWcsException
     * @throws URISyntaxException
     * @throws FileNotFoundException
     */
    public static JWcsMap readFits(String url, int hdu) throws JWcsException, URISyntaxException, 
        FileNotFoundException {
        final Map<String, String> keyMap = new HashMap<>();
        final URI uri = new URI(url);
        try (Fits fits = new Fits(uri.toURL())) {
            final Header hdr = fits.getHDU(hdu).getHeader();
            final Cursor<String, HeaderCard> c = hdr.iterator();
            while (c.hasNext()) {
                final HeaderCard card = c.next();
                keyMap.put(card.getKey(), card.getValue());
            }
        } catch (NoSuchMethodError | Exception ex) {
            ex.printStackTrace();
            final HeaderFitsReader hdr = new HeaderFitsReader(uri);
            final List<List<String>> listKeywords = hdr.readKeywords();
            listKeywords.stream().forEach(keywordLine -> 
                keyMap.put(keywordLine.get(0), keywordLine.get(1))
            );
        }
        final JWcsMap wcs = new JWcsMap(keyMap);
        wcs.doInit();
        return wcs;
    }
    
    /**
     * Reads a fits file and returns the {@linkplain JWcsMap} object to perform 
     * pixel <--> sky coordinate conversions.
     * @param file The fits file to read.
     * @param hdu HDU number.
     * @return The instance of {@linkplain JWcsMap}.
     * @throws JWcsException
     * @throws FileNotFoundException
     */
    public static JWcsMap readFits(File file, int hdu) throws JWcsException, 
        FileNotFoundException {
        final Map<String, String> keyMap = new HashMap<>();
        try (Fits fits = new Fits(file)) {
            final Header hdr = fits.getHDU(hdu).getHeader();
            final Cursor<String, HeaderCard> c = hdr.iterator();
            while (c.hasNext()) {
                final HeaderCard card = c.next();
                keyMap.put(card.getKey(), card.getValue());
            }
        } catch (NoSuchMethodError | Exception ex) {
            ex.printStackTrace();
            final HeaderFitsReader hdr = new HeaderFitsReader(file);
            final List<List<String>> listKeywords = hdr.readKeywords();
            listKeywords.stream().forEach(keywordLine -> 
                keyMap.put(keywordLine.get(0), keywordLine.get(1))
            );
        }
        final JWcsMap wcs = new JWcsMap(keyMap);
        wcs.doInit();
        return wcs;
    }

    /**
     * Reads a String with degrees, minutes, and arcseconds to get those 
     * individual fields. Format can be ddd mm ss.ss, ddd:mm:ss.ss, or d.ddd.
     * @param dms The String with degrees, minutes, and arcseconds.
     * @return The instance of {@linkplain DMS}.
     */
    public static DMS getDMS(String dms) {
        return new DMS(dms);
    }
    
    /**
     * Reads a String with hours, minutes, and seconds of time to get those 
     * individual fields. Format can be hh mm ss.ss, hh:mm:ss.ss, hh, or d.ddd.
     * IMPORTANT: In case of d.ddd input, with decimal point, the value is 
     * assumed to be in degrees. Otherwise, for a integer value, is assumed 
     * to be in hours.
     * @param dms The String with degrees, minutes, and arcseconds.
     * @return The instance of {@linkplain HMS}.
     */
    public static HMS getHMS(String hms) {
        return new HMS(hms);
    }
}
