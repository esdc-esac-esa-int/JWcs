/*
 * Copyright (C) 2016 malapert
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package io.github.malapert.jwcs.proj;

import io.github.malapert.jwcs.proj.exception.BadProjectionParameterException;
import io.github.malapert.jwcs.utility.NumericalUtility;
import java.util.logging.Logger;
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
        AbstractZenithalProjection instance = proj;
        double expResult = 0.0;
        double result = instance.getPhi0();
        assertEquals(expResult, result, 0.0);
    }

    /**
     * Test of getTheta0 method, of class AbstractZenithalProjection.
     */
    @Test
    public void testGetTheta0() {
        System.out.println("getTheta0");
        AbstractZenithalProjection instance = proj;
        double expResult = NumericalUtility.HALF_PI;
        double result = instance.getTheta0();
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
        double expResult = phi0;
        double result = instance.getPhi0();        
        assertEquals(expResult, result, 0.001);
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
        AbstractZenithalProjection instance = proj;
        String expResult = "Zenithal (azimuthal) projections";
        String result = instance.getNameFamily();
        assertEquals(expResult, result);
    }  
}
