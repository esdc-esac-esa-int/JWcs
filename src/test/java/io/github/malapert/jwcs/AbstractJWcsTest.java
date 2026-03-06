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

import io.github.malapert.jwcs.crs.Equatorial;
import io.github.malapert.jwcs.crs.AbstractCrs;
import io.github.malapert.jwcs.proj.exception.JWcsException;
import java.io.IOException;
import java.net.URL;
import nom.tam.fits.Fits;
import nom.tam.fits.FitsException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Jean-Christophe Malapert
 */
public class AbstractJWcsTest {
    
    AbstractJWcs wcs;
    AbstractJWcs wcs1;
    
    /**
     *
     * @throws FitsException
     * @throws IOException
     * @throws JWcsException
     */
    public AbstractJWcsTest() throws FitsException, IOException, JWcsException {

        URL url = getClass().getClassLoader().getResource("WFPC2ASSNu5780205bx.fits");
        assertNotNull(url);
        Fits simpleFits = new Fits(url.toString());
        wcs = new JWcsFits(simpleFits);
        wcs.doInit();
        url = getClass().getClassLoader().getResource("FOCx38i0101t_c0f.fits");
        assertNotNull(url);
        Fits fits = new Fits(url.toString());
        wcs1 = new JWcsFits(fits, 0);
        wcs1.doInit();
    }
    
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


    /**
     * Test of getCrs method, of class AbstractJWcs.
     */
    @Test
    public void testGetCrs() {
        System.out.println("getSkySystem");
        AbstractCrs result = wcs.getCrs();
        double equinox = ((Equatorial)result).getEquinox();
        Double epoch = ((Equatorial)result).getEpochObs();
        String refSystem = ((Equatorial)result).getCoordinateReferenceFrame().getReferenceFrame().name();       
        assertEquals(2000.0, equinox, 1e-12);
        if (Double.isNaN(epoch)) {
            epoch = null;
        }
        assertNull(epoch);
        assertEquals("FK5", refSystem);
        
        result = wcs1.getCrs();
        equinox = ((Equatorial)result).getEquinox();
        epoch = ((Equatorial)result).getEpochObs();
        refSystem = result.getCoordinateReferenceFrame().getReferenceFrame().name();
        assertEquals(2000.0, equinox, 1e-12);
        if (Double.isNaN(epoch)) {
            epoch = null;
        }
        assertNull(epoch);
        assertEquals("FK5", refSystem);        
    }

    /**
     * Test of pix2wcs method, of class AbstractJWcs.
     * @throws java.lang.Exception
     */
    @Test
    public void testPix2wcsDoubleArr() throws Exception {
        System.out.println("pix2wcs");
        double[] pixels = new double[]{1,1,1,100,100,100,100,1};
        double[] expResult = new double[]{182.632805,39.406235,182.632538,39.403504,182.636069,39.403298,182.636336,39.406029};
        double[] result = wcs.pix2wcs(pixels);
        assertArrayEquals(expResult, result, 1e-6);
        
        pixels = new double[]{1,1,1,1024,1024,1024,1024,1};
        expResult = new double[]{182.64307,39.406175,182.63451,39.411264,182.62843,39.405157,182.63699,39.40007};
        result = wcs1.pix2wcs(pixels);
        assertArrayEquals(expResult, result, 1e-5);        
    }

    /**
     * Test of getCenter method, of class AbstractJWcs.
     * @throws java.lang.Exception
     */
    @Test
    public void testGetCenter() throws Exception {
        System.out.println("getCenter");
        final double[] expResult = new double[]{182.63442,39.404782};
        final double[] result = wcs.getCenter();
        assertArrayEquals(expResult, result, 1e-5);
    }

    /**
     * Test of getFov method, of class AbstractJWcs.
     * @throws java.lang.Exception
     */
    @Test
    public void testGetFov() throws Exception {
        System.out.println("getFov");
        final double[] pixels = new double[]{0.5,0.5,100+0.5,0.5,100+0.5,100+0.5,0.5,100+0.5};
        final double[] expResult = wcs.pix2wcs(pixels);
        final double[] result = wcs.getFov();
        assertArrayEquals(expResult, result, 1e-5);
    }    
}
