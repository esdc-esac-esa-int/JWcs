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

import io.github.malapert.jwcs.JWcsFits;
import io.github.malapert.jwcs.proj.exception.JWcsException;
import io.github.malapert.jwcs.proj.exception.ProjectionException;
import nom.tam.fits.Fits;
import nom.tam.fits.FitsException;
import org.junit.*;

import java.io.IOException;
import java.util.Objects;

import static org.junit.Assert.assertArrayEquals;

/**
 * ZPN unit test.
 *
 * @author Jean-Christophe Malapert
 */
public class ZPNTest extends AbstractProjectionTest {

    /**
     *
     * @throws FitsException
     * @throws IOException
     * @throws JWcsException
     */
    public ZPNTest() throws FitsException, IOException, JWcsException {
        super(new JWcsFits(new Fits(
                Objects.requireNonNull(ZPNTest.class.getClassLoader().getResource("1904-66_ZPN.fits")).toString())));
    }

    @BeforeClass
    public static void setUpClass() {
        // do nothing        
    }

    @AfterClass
    public static void tearDownClass() {
        // do nothing
    }

    @Before
    public void setUp() {
        // do nothing

    }

    @After
    public void tearDown() {
        // do nothing
    }

    /**
     * Test of project method, of class ZPN.
     *
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException
     */
    @Test
    public void testProjectZPN() throws ProjectionException {
        System.out.println("project ZPN");
        final double[][] expectedResults =
                { { 263.471000708007352, -78.497682328997385 }, { 266.783268968517177, -50.24524022078576 },
                  { 294.357836271455028, -39.77023899472649 }, { 312.674220190438405, -71.468154470727242 } };
        double[] result = wcs.pix2wcs(1, 1);
        assertArrayEquals(expectedResults[0], result, 1e-13);

        result = wcs.pix2wcs(192, 1);
        assertArrayEquals(expectedResults[1], result, 1e-13);

        result = wcs.pix2wcs(192, 192);
        assertArrayEquals(expectedResults[2], result, 1e-13);

        result = wcs.pix2wcs(1, 192);
        assertArrayEquals(expectedResults[3], result, 1e-12);
    }

    /**
     * Test of projectInverse method, of class ZPN.
     *
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException
     */
    @Test
    public void testProjectInverseZPN() throws ProjectionException {
        System.out.println("projectInverse ZPN");
        final double[][] expectedResults = { { 1.0d, 1.0d }, { 192.d, 1.0d }, { 192.d, 192d }, { 1.0d, 192d } };
        double[] result;
        for (final double[] expectedResult : expectedResults) {
            result = wcs.pix2wcs(expectedResult);
            result = wcs.wcs2pix(result);
            assertArrayEquals(expectedResult, result, 1e-12);
        }
    }
}
