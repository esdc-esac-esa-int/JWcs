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
import org.junit.Test;

import java.io.IOException;
import java.util.Objects;

import static org.junit.Assert.assertArrayEquals;

/**
 * CEA unit test.
 * @author Jean-Christophe Malapert
 */
public class CEATest extends AbstractProjectionTest {

    /**
     *
     * @throws JWcsException
     * @throws FitsException
     * @throws IOException
     */
    public CEATest() throws JWcsException, FitsException, IOException {
        super(new JWcsFits(new Fits(
                Objects.requireNonNull(CEATest.class.getClassLoader().getResource("1904-66_CEA.fits")).toString())));
    }

    /**
     *
     * @throws ProjectionException
     */
    @Test
    public void testProjectCEA() throws ProjectionException {
        System.out.println("project CEA");
        final double[][] expectedResults =
                { { 268.440852654621153, -73.379693805485672 }, { 269.090243859417001, -60.649088748117173 },
                  { 294.131910549115673, -58.362095662786679 }, { 307.520448192392109, -69.383029011108476 } };
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
     *
     * @throws ProjectionException
     */
    @Test
    public void testProjectInverseCEA() throws ProjectionException {
        System.out.println("projectInverse CEA");
        final double[][] expectedResults = { { 1.0d, 1.0d }, { 192.d, 1.0d }, { 192.d, 192d }, { 1.0d, 192d } };
        double[] result;
        for (final double[] expectedResult : expectedResults) {
            result = wcs.pix2wcs(expectedResult);
            result = wcs.wcs2pix(result);
            assertArrayEquals(expectedResult, result, 1e-12);
        }
    }

}
