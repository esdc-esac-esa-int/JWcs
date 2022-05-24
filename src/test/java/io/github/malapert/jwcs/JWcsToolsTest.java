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

import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import io.github.malapert.jwcs.JWcsTools.ANGLE;
import io.github.malapert.jwcs.JWcsTools.COORD_SYS;
import io.github.malapert.jwcs.JWcsTools.FRAME;
import io.github.malapert.jwcs.position.SkyPosition;
import io.github.malapert.jwcs.utility.DMS;
import io.github.malapert.jwcs.utility.HMS;

/**
 *
 * @author T. Alonso Albi (tomas.alonsoalbi@esa.int)
 */
public class JWcsToolsTest {

    /**
    *
    */
   @BeforeClass
   public static void setUpClass() {
       //do nothing
   }
   
   /**
    *
    */
   @AfterClass
   public static void tearDownClass() {
       //do nothing
   }
   
   /**
    *
    */
   @Before
   public void setUp() {
       //do nothing
   }
   
   /**
    *
    */
   @After
   public void tearDown() {
       //do nothing
   }
   
   @Test
   public void convertAngleTest() {
       assertEquals(1, ANGLE.DEG.getConversionFactorTo(ANGLE.DEG), 1E-12);
       assertEquals(60, ANGLE.DEG.getConversionFactorTo(ANGLE.ARCMIN), 1E-12);
       assertEquals(3600, ANGLE.DEG.getConversionFactorTo(ANGLE.ARCSEC), 1E-12);
       assertEquals(1.0 / 15.0, ANGLE.DEG.getConversionFactorTo(ANGLE.HOUR), 1E-12);
   }
   
   @Test
   public void convertSkyPositionTest() {
       SkyPosition pos = JWcsTools.getSkyPosition(90, 0, COORD_SYS.EQUATORIAL);
       SkyPosition out = JWcsTools.convert(pos, COORD_SYS.ECLIPTIC);
       assertEquals(90, out.getLongitude(), 1E-12);
       assertEquals(-23.4392911111111, out.getLatitude(), 1E-12);
   }

   @Test
   public void skySeparationTest() {
       SkyPosition pos1 = JWcsTools.getSkyPosition(0, 0, COORD_SYS.EQUATORIAL);
       SkyPosition pos2 = JWcsTools.getSkyPosition(90, 0, COORD_SYS.EQUATORIAL);
       assertEquals(90, JWcsTools.skySeparation(pos1, pos2), 1E-12);
   }
   
   @Test
   public void precessTest() {
       double ra = 0;
       double dec = 0;
       double out[] = JWcsTools.precessEquinox(ra, dec, 2000, 2050, FRAME.FK5);
       assertEquals(0.6407181591449038, out[0], 1E-12);
       assertEquals(0.2783410827199949, out[1], 1E-12);
   }
   
   @Test
   public void convertB1950Test() {
       double out1[] = JWcsTools.convertB1950ToICRS(0, 0);
       double out2[] = JWcsTools.convertB1950ToFK5(0, 0);
       assertEquals(0.6406846124680589, out1[0], 1E-12);
       assertEquals(0.27840696926887554, out1[1], 1E-12);
       assertEquals(0.6406910005754157, out2[0], 1E-12);
       assertEquals(0.2784094350773699, out2[1], 1E-12);
   }
   
   @Test
   public void convertJDAndISOdatesTest() throws Exception {
       double jd = 2456915;
       String iso = "2014-09-14T12:00:00";
       
       assertEquals(jd, JWcsTools.isoToJulianDate(iso), 1E-12);
       assertEquals(iso, JWcsTools.julianDateToISO(jd));
   }
   
   @Test
   public void convertDMSandHMSTest() throws Exception {
       String sdms = "10 20 30";
       String shms = "-15.2";
       DMS dms = JWcsTools.getDMS(sdms);
       HMS hms = JWcsTools.getHMS(shms);
       
       assertEquals(10, dms.getDegrees(), 1E-12);
       assertEquals(20, dms.getMin(), 1E-12);
       assertEquals(30, dms.getSec(), 1E-12);
       assertEquals(1, hms.getHours(), 1E-12);
       assertEquals(0, hms.getMin(), 1E-12);
       assertEquals(48, hms.getSec(), 1E-6);
       assertEquals(-1, hms.getSign(), 1E-12);
   }
   
   @Test
   public void readFitsTest() throws Exception {
       JWcsMap wcs = JWcsTools.readFits("http://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.fits", 0);
       double sky[] = wcs.pix2wcs(0, 0);
       assertEquals(52.53373984070186, sky[0], 1E-12);
       assertEquals(-28.760675854311447, sky[1], 1E-12);
   }
}
