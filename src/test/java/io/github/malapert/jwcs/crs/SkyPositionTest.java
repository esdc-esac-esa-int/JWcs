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

import io.github.malapert.jwcs.datum.FK5;
import io.github.malapert.jwcs.datum.ICRS;
import io.github.malapert.jwcs.position.SkyPosition;
import org.junit.*;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 *
 * @author malapert
 */
public class SkyPositionTest {

    public SkyPositionTest() {
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
     * Test of getLongitude method, of class SkyPosition.
     */
    @Test
    public void testGetLongitude() {
        System.out.println("getLongitude");
        SkyPosition instance = new SkyPosition(10, 0, new Equatorial());
        double expResult = 10.0;
        double result = instance.getLongitude();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of getLongitudeAsSexagesimal method, of class SkyPosition.
     */
    @Test
    public void testGetLongitudeAsSexagesimal() {
        System.out.println("getLongitudeAsSexagesimal");
        SkyPosition instance = new SkyPosition(10.5, 0, new Equatorial());
        String expResult = "00:42:00.000";
        String result = instance.getLongitudeAsSexagesimal();
        assertEquals(expResult, result);
    }

    /**
     * Test of getLatitude method, of class SkyPosition.
     */
    @Test
    public void testGetLatitude() {
        System.out.println("getLatitude");
        SkyPosition instance = new SkyPosition(10.5, 0, new Equatorial());
        double expResult = 0.0;
        double result = instance.getLatitude();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of getLatitudeAsSexagesimal method, of class SkyPosition.
     */
    @Test
    public void testGetLatitudeAsSexagesimal() {
        System.out.println("getLatitudeAsSexagesimal");
        SkyPosition instance = new SkyPosition(10.5, 0, new Equatorial());
        String expResult = "+00:00:00.00";
        String result = instance.getLatitudeAsSexagesimal();
        assertEquals(expResult, result);
    }

    /**
     * Test of setLongitude method, of class SkyPosition.
     */
    @Test
    public void testSetLongitude() {
        System.out.println("setLongitude");
        double longitude = 20.0;
        SkyPosition instance = new SkyPosition(0, 0, new Equatorial());
        instance.setLongitude(longitude);
        double expResult = 20.0;
        double result = instance.getLongitude();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of setLatitude method, of class SkyPosition.
     */
    @Test
    public void testSetLatitude() {
        System.out.println("setLatitude");
        double latitude = -20.0;
        SkyPosition instance = new SkyPosition(0, 0, new Equatorial());
        instance.setLatitude(latitude);
        double expResult = -20.0;
        double result = instance.getLatitude();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of getDoubleArray method, of class SkyPosition.
     */
    @Test
    public void testGetDoubleArray() {
        System.out.println("getDoubleArray");
        SkyPosition instance = new SkyPosition(10, -20, new Equatorial());
        double[] expResult = new double[] { 10, -20 };
        double[] result = instance.getDoubleArray();
        assertArrayEquals(expResult, result, 1e-13);
    }

    /**
     * Test of equals method, of class SkyPosition.
     */
    @Test
    public void testEquals() {
        System.out.println("equals");
        Object obj = new SkyPosition(0, 0, new Equatorial(new ICRS()));
        Object obj2 = new SkyPosition(0, 0, new Equatorial(new FK5()));
        SkyPosition instance = new SkyPosition(0, 0, new Equatorial());
        boolean expResult = true;
        boolean expResult2 = false;
        boolean result = instance.equals(obj);
        boolean result2 = instance.equals(obj2);
        assertEquals(expResult, result);
        assertEquals(expResult2, result2);

    }
}
