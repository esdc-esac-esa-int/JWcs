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

import io.github.malapert.jwcs.datum.ICRS;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author malapert
 */
public class ICRSTest {
    
    public ICRSTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }


    /**
     * Test of getEpochObs method, of class ICRS.
     */
    @Test
    public void testGetEpochObs() {
        System.out.println("getEpochObs");
        ICRS instance = new ICRS();
        double expResult = Double.NaN;
        double result = instance.getEpochObs();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of getEquinox method, of class ICRS.
     */
    @Test
    public void testGetEquinox() {
        System.out.println("getEquinox");
        ICRS instance = new ICRS();
        double expResult = 2000.0;
        double result = instance.getEquinox();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of setEquinox method, of class ICRS.
     */
    @Test
    public void testSetEquinox_String() {
        System.out.println("setEquinox");
        String equinox = "J2001";
        ICRS instance = new ICRS();
        instance.setEquinox(equinox);
        double expResult = 2000.0;
        double result = instance.getEquinox();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of setEquinox method, of class ICRS.
     */
    @Test
    public void testSetEquinox_double() {
        System.out.println("setEquinox");
        double equinox = 1950;
        ICRS instance = new ICRS();
        instance.setEquinox(equinox);
        double expResult = 2000;
        double result = instance.getEquinox();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of setEpochObs method, of class ICRS.
     */
    @Test
    public void testSetEpochObs_String() {
        System.out.println("setEpochObs");
        String epochObs = "B1950";
        ICRS instance = new ICRS();
        instance.setEpochObs(epochObs);
        double expResult = Double.NaN;
        double result = instance.getEpochObs();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of setEpochObs method, of class ICRS.
     */
    @Test
    public void testSetEpochObs_double() {
        System.out.println("setEpochObs");
        double epochObs = 1950.0;
        ICRS instance = new ICRS();
        instance.setEpochObs(epochObs);
        double expResult = Double.NaN;
        double result = instance.getEpochObs();
        assertEquals(expResult, result, 0.0);
    }


    /**
     * Test of equals method, of class ICRS.
     */
    @Test
    public void testEquals() {
        System.out.println("equals");
        ICRS instance1 = new ICRS();
        ICRS instance2 = new ICRS();
        boolean result = instance1.equals(instance2);
        assertTrue(result);
    }
    
}
