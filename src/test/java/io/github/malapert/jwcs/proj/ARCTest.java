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
 * ARC unit test.
 * @author Jean-Christophe Malapert
 */
public class ARCTest extends AbstractProjectionTest {
    
    /**
     *
     * @throws FitsException
     * @throws IOException
     * @throws JWcsException
     */
    public ARCTest() throws FitsException, IOException, JWcsException {
        super(new JWcsFits(new Fits(new URL("http://tdc-www.harvard.edu/wcstools/samples/1904-66_ARC.fits"))));
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
     * Test of project method, of class ARC.
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException
     */
    @Test
    public void testProjectARC() throws ProjectionException {
        System.out.println("project ARC for particular points");
        final double expectedResults[][] = {
            {269.056730777738039,  -73.468299585347012},
            {269.467149632953806,  -60.735941026372636},
            {293.066101937638564,  -58.194463838114913},
            {307.011804331818496,  -69.29965938606621 }
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
     * Test of projectInverse method, of class AIT.
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException
     */
    @Test
    public void testProjectInverseARC() throws ProjectionException {
        System.out.println("projectInverse ARC");
        final double expectedResults[][] = {
            {1.0d, 1.0d},
            {192.d, 1.0d},
            {192.d, 192d},
            {1.0d, 192d}
        };   
        double[] result;
        for (final double[] expectedResult : expectedResults) {
            result = wcs.pix2wcs(expectedResult);
            result = wcs.wcs2pix(result);
             assertArrayEquals(expectedResult, result, 1e-12);
        }        
    }
    
}
