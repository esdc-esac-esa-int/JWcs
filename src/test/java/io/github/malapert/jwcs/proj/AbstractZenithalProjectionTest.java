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
package io.github.malapert.jwcs.proj;

import io.github.malapert.jwcs.proj.exception.BadProjectionParameterException;
import io.github.malapert.jwcs.utility.NumericalUtility;
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
public class AbstractZenithalProjectionTest {
    
    private final AbstractZenithalProjection proj;
    
    public AbstractZenithalProjectionTest() throws BadProjectionParameterException {
        this.proj = new AZP();
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
     * Test of getPhi0 method, of class AbstractZenithalProjection.
     */
    @Test
    public void testGetPhi0() {
        System.out.println("getPhi0");
        double expResult = 0.0;
        double result = proj.getPhi0();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of getTheta0 method, of class AbstractZenithalProjection.
     */
    @Test
    public void testGetTheta0() {
        System.out.println("getTheta0");
        double expResult = NumericalUtility.HALF_PI;
        double result = proj.getTheta0();
        assertEquals(expResult, result, 0.001);
    }

    /**
     * Test of setPhi0 method, of class AbstractZenithalProjection.
     */
    @Test
    public void testSetPhi0() {
        System.out.println("setPhi0");
        double phi0 = 2.0;
        AbstractZenithalProjection instance = proj;
        instance.setPhi0(phi0);
        double result = instance.getPhi0();
        assertEquals(phi0, result, 0.001);
    }

    /**
     * Test of setTheta0 method, of class AbstractZenithalProjection.
     */
    @Test
    public void testSetTheta0() {
        System.out.println("setTheta0");
        double theta0 = 4.0;
        AbstractZenithalProjection instance = proj;
        instance.setTheta0(theta0);
        double expResult = 4.0;
        double result = instance.getTheta0();
        assertEquals(expResult, result, 0.001);        
    }

    /**
     * Test of getNameFamily method, of class AbstractZenithalProjection.
     */
    @Test
    public void testGetNameFamily() {
        System.out.println("getNameFamily");
        String expResult = "Zenithal (azimuthal) projections";
        String result = proj.getNameFamily();
        assertEquals(expResult, result);
    }  
}
