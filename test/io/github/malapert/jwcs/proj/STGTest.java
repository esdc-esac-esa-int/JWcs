/* 
 * Copyright (C) 2014-2016 Jean-Christophe Malapert
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

import io.github.malapert.jwcs.proj.exception.JWcsException;
import io.github.malapert.jwcs.JWcsFits;
import io.github.malapert.jwcs.proj.exception.ProjectionException;
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
 * STG unit test.
 * @author Jean-Christophe Malapert
 */
public class STGTest extends AbstractProjectionTest {
    
    /**
     *
     * @throws FitsException
     * @throws IOException
     * @throws JWcsException
     */
    public STGTest() throws FitsException, IOException, JWcsException {
        super(new JWcsFits(new Fits(new URL("http://tdc-www.harvard.edu/wcstools/samples/1904-66_STG.fits"))));
    }
    
    /**
     *
     */
    @BeforeClass
    public static void setUpClass() {
    }
    
    /**
     *
     */
    @AfterClass
    public static void tearDownClass() {
    }
    
    /**
     *
     */
    @Before
    public void setUp() {
    }
    
    /**
     *
     */
    @After
    public void tearDown() {
    }

    /**
     * Test of project method, of class STG.
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException
     */
    @Test
    public void testProjectSTG() throws ProjectionException {
        System.out.println("project STG");
        double expectedResults[][] = {
            { 269.37825680266144 ,  -73.256130460250546},
            { 269.645741550822208,  -61.036020984514579},
            { 292.979346455154086,  -58.65820590406657 },
            { 306.658467125869947,  -69.210340276541672}
        };
        double[] result = wcs.pix2wcs(1, 1);
        assertArrayEquals(expectedResults[0], result, 1e-13);

        result = wcs.pix2wcs(192, 1);
        assertArrayEquals(expectedResults[1], result, 1e-13);

        result = wcs.pix2wcs(192, 192);
        assertArrayEquals(expectedResults[2], result, 1e-13);

        result = wcs.pix2wcs(1, 192);
        assertArrayEquals(expectedResults[3], result, 1e-13);
    }

    /**
     * Test of projectInverse method, of class STG.
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException
     */
    @Test
    public void testProjectInverseSTG() throws ProjectionException {
        System.out.println("projectInverse STG");
        double expectedResults[][] = {
            {1.0d, 1.0d},
            {192.d, 1.0d},
            {192.d, 192d},
            {1.0d, 192d}
        };   
        double[] result;
        for (double[] expectedResult : expectedResults) {
            result = wcs.pix2wcs(expectedResult);
            result = wcs.wcs2pix(result);
             assertArrayEquals(expectedResult, result, 1e-12);
        }  
    }
    
}
