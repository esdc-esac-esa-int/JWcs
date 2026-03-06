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
package io.github.malapert.jwcs.utility;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Tests of Utility package.
 * @author Jean-Christophe Malapert
 */
public class TimeUtilityTest {
    
    public TimeUtilityTest() {
        //do nothing
    }
    
    @BeforeClass
    public static void setUpClass() {
        //do nothing
    }
    
    @AfterClass
    public static void tearDownClass() {
        //do nothing
    }
    
    @Before
    public void setUp() {
        //do nothing
    }
    
    @After
    public void tearDown() {
        //do nothing
    }


    /**
     * Test of convertJulianDateToISO method, of class TimeUtility.
     */
    @Test
    public void testConvertJulianDateToISO() {
        System.out.println("julianDateToISO");
        final double julianDate = 2456915;
        final String expResult = "2014-09-14T12:00:00";
        final String result = TimeUtility.convertJulianDateToISO(julianDate);
        assertEquals(expResult, result);
    }

    /**
     * Test of convertModifiedJulianDateToISO method, of class TimeUtility.
     */
    @Test
    public void testConvertModifiedJulianDateToISO() {
        System.out.println("modifiedJulianDateToISO");
        final double modifiedJulianDate = 53005;
        final String expResult = "2004-01-01T00:00:00";
        final String result = TimeUtility.convertModifiedJulianDateToISO(modifiedJulianDate);
        assertEquals(expResult, result);
    }

    /**
     * Test of convertISOToJulianDate method, of class TimeUtility.
     * @throws java.lang.Exception
     */
    @Test
    public void testConvertISOToJulianDate() throws Exception {
        System.out.println("ISOToJulianDate");
        final String dateObs = "2014-09-14T12:08:18";
        final double expResult = 2456915.005763889;
        final double result = TimeUtility.convertISOToJulianDate(dateObs);
        assertEquals(expResult, result, 1e-9);             
    }
    
    
    /**
     * Test of convertISOToModifiedJulianDate method, of class TimeUtility.
     * @throws java.lang.Exception
     */
    @Test
    public void testConvertISOToModifiedJulianDate() throws Exception {
        System.out.println("ISOToModifiedJulianDate");
        
        final String datObs = "2004-01-01T00:00:00.000";        
        final double result = TimeUtility.convertISOToModifiedJulianDate(datObs);
        final double expResultMJD = 53005.0;
        assertEquals(expResultMJD, result, 1e-6);
    }
    
    /**
     *
     * @throws Exception
     */
    @Test
    public void convertAndreverse() throws Exception {
        System.out.println("Convert and reverse");
        
        String dateObs = "10/03/78";        
        double result = TimeUtility.convertISOToModifiedJulianDate(dateObs);
        String expResult = TimeUtility.convertModifiedJulianDateToISO(result);
        assertEquals("1978-03-10T00:00:00", expResult);
        
        dateObs = "1978-03-10";        
        result = TimeUtility.convertISOToModifiedJulianDate(dateObs);
        expResult = TimeUtility.convertModifiedJulianDateToISO(result);
        assertEquals("1978-03-10T00:00:00", expResult);
        
        dateObs = "1978-03-10T10:20:02";        
        result = TimeUtility.convertISOToModifiedJulianDate(dateObs);
        expResult = TimeUtility.convertModifiedJulianDateToISO(result);
        assertEquals(dateObs, expResult);
    }

    /**
     * Test of convertEpochJulian2JD method, of class TimeUtility.
     */
    @Test
    public void testConvertEpochJulian2JD() {
        System.out.println("epochJulian2JD");
        final double jEpoch = 1983.99863107d;
        final double expResult = 2445700.5d;
        final double result = TimeUtility.convertEpochJulian2JD(jEpoch);
        assertEquals(expResult, result, 1e-1);
    }

    /**
     * Test of convertJD2epochJulian method, of class TimeUtility.
     */
    @Test
    public void testConvertJD2epochJulian() {
        System.out.println("JD2epochJulian");
        final double jd = 2445700.5d;
        final double expResult = 1983.99863107d;
        final double result = TimeUtility.convertJD2epochJulian(jd);
        assertEquals(expResult, result, 1e-7);
    }

    /**
     * Test of convertEpochBessel2JD method, of class TimeUtility.
     */
    @Test
    public void testConvertEpochBessel2JD() {
        System.out.println("epochBessel2JD");
        final double bEpoch = 1983.99956681d;
        final double expResult = 1983.99956681d;
        final double jd = TimeUtility.convertEpochBessel2JD(bEpoch);
        final double result = TimeUtility.convertJD2epochBessel(jd);
        assertEquals(expResult, result, 1e-7);
    }

    /**
     * Test of epochs method, of class TimeUtility.
     */
    @Test
    public void testEpochs() {
        System.out.println("epochs");
        final String epoch = "F2008-03-31T8:09";
        final double[] expResult = new double[]{2008.2474210134737, 2008.2459673739454, 2454556.8395833336};
        final double[] result = TimeUtility.epochs(epoch);
        assertArrayEquals(expResult, result, 1e-12);
    }

    /**
     * Test of epochs method, of class TimeUtility.
     */
    @Test
    public void testEpochs1() {
        System.out.println("epochs");
        final String epoch = "F2007-01-14T13:18:59.9";
        final double[] expResult = new double[]{2007.0378545262108, 2007.0364267212976, 2454115.0548599539};
        final double[] result = TimeUtility.epochs(epoch);
        assertArrayEquals(expResult, result, 1e-12);
    }

    /**
     * Test of epochs method, of class TimeUtility.
     */
    @Test
    public void testEpochs2() {
        System.out.println("epochs");
        final String epoch = "j2007.0364267212976";
        final double[] expResult = new double[]{2007.0378545262108, 2007.0364267212976, 2454115.0548599539};
        final double[] result = TimeUtility.epochs(epoch);
        assertArrayEquals(expResult, result, 1e-12);
    }

    /**
     * Test of epochs method, of class TimeUtility.
     */
    @Test
    public void testEpochs3() {
        System.out.println("epochs");
        final String epoch = "b2007.0378545262108";
        final double[] expResult = new double[]{2007.0378545262108, 2007.0364267212976, 2454115.0548599539};
        final double[] result = TimeUtility.epochs(epoch);
        assertArrayEquals(expResult, result, 1e-12);
    }
    

    /**
     * Test of jd method, of class TimeUtility.
     */
    @Test
    public void testJd() {
        System.out.println("jd");
        final int year = -4712;
        final int month = 1;
        final double dayNumber = 1.5d;
        final double expResult = 0.0d;
        final double result = TimeUtility.jd(year, month, dayNumber);
        assertEquals(expResult, result, 0.0);
    }
    
    /**
     * Test of jd method, of class TimeUtility.
     */
    @Test
    public void testJd1() {
        System.out.println("jd");
        final int year = 0;
        final int month = 1;
        final double dayNumber = 1.0d;
        final double expResult = 1721057.5d;
        final double result = TimeUtility.jd(year, month, dayNumber);
        assertEquals(expResult, result, 0.0);
    }    

    /**
     * Test of jd method, of class TimeUtility.
     */
    @Test
    public void testJd2() {
        System.out.println("jd");
        final int year = 1582;
        final int month = 10;
        final double dayNumber = 4.0d;
        final double expResult = 2299159.5d;
        final double result = TimeUtility.jd(year, month, dayNumber);
        assertEquals(expResult, result, 0.0);
    }      
    
    /**
     * Test of jd method, of class TimeUtility.
     */
    @Test
    public void testJd3() {
        System.out.println("jd");
        final int year = 1582;
        final int month = 10;
        final double dayNumber = 15.0d;
        final double expResult = 2299170.5d;
        final double result = TimeUtility.jd(year, month, dayNumber);
        assertEquals(expResult, result, 0.0);
    }       
    
    /**
     * Test of jd method, of class TimeUtility.
     */
    @Test
    public void testJd4() {
        System.out.println("jd");
        final int year = 1582;
        final int month = 10;
        final double dayNumber = 15.5d;
        final double expResult = 2299161.0d;
        final double result = TimeUtility.jd(year, month, dayNumber);
        assertEquals(expResult, result, 0.0);
    }     
    
    /**
     * Test of jd method, of class TimeUtility.
     */
    @Test
    public void testJd5() {
        System.out.println("jd");
        final int year = 1970;
        final int month = 1;
        final double dayNumber = 1d;
        final double expResult = 2440587.5d;
        final double result = TimeUtility.jd(year, month, dayNumber);
        assertEquals(expResult, result, 0.0);
    }      
}
