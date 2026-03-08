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
 * SFL unit test.
 * @author malapert
 */
public class SFLTest extends AbstractProjectionTest {

    /**
     *
     * @throws FitsException
     * @throws IOException
     * @throws JWcsException
     */
    public SFLTest() throws FitsException, IOException, JWcsException {
        super(new JWcsFits(new Fits(
                Objects.requireNonNull(SFLTest.class.getClassLoader().getResource("1904-66_SFL.fits")).toString())));
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
     * Test of project method, of class SFL.
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException
     */
    @Test
    public void testProjectSFL() throws ProjectionException {
        System.out.println("project SFL");
        final double[][] expectedResults =
                { { 268.467379871114304, -73.504056521463767 }, { 269.108796034278896, -60.772983212159005 },
                  { 293.614959948683349, -57.878452615484491 }, { 306.902910888361077, -69.223289162578311 } };
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
     * Test of projectInverse method, of class SFL.
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException
     */
    @Test
    public void testProjectInverseSFL() throws ProjectionException {
        System.out.println("projectInverse SFL");
        final double[][] expectedResults = { { 1.0d, 1.0d }, { 192.d, 1.0d }, { 192.d, 192d }, { 1.0d, 192d } };
        double[] result;
        for (final double[] expectedResult : expectedResults) {
            result = wcs.pix2wcs(expectedResult);
            result = wcs.wcs2pix(result);
            assertArrayEquals(expectedResult, result, 1e-12);
        }
    }

}
