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

import io.github.malapert.jwcs.position.SkyPosition;
import io.github.malapert.jwcs.proj.exception.JWcsException;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.util.logging.Level;
import java.util.logging.Logger;

import static org.junit.Assert.assertEquals;

/**
 * Projection test.
 *
 * @author Jean-Christophe Malapert
 */
@Ignore
public class ConverterTest {

    /**
     *
     */
    protected AbstractCrs source;

    /**
     *
     */
    protected AbstractCrs target;

    /**
     *
     */
    protected double minLat;

    /**
     *
     */
    protected double maxLat;

    /**
     *
     */
    protected double tolerance;

    private static final double TOLERANCE = 1.0e-7;
    static final Logger LOG = Logger.getLogger("");

    /**
     *
     */
    public ConverterTest() {
        LOG.setLevel(Level.OFF);
        this.tolerance = TOLERANCE;
    }

    /**
     *
     * @param source
     * @param target
     * @throws JWcsException
     */
    public ConverterTest(final AbstractCrs source, final AbstractCrs target) throws JWcsException {
        this(source, target, TOLERANCE);
    }

    /**
     *
     * @param source
     * @param target
     * @param tolerance
     */
    public ConverterTest(final AbstractCrs source, final AbstractCrs target, final double tolerance) {
        this.source = source;
        this.target = target;
        this.tolerance = tolerance;
        LOG.setLevel(Level.OFF);
    }

    /**
     *
     */
    @BeforeClass
    public static void setUpClass() {
        System.out.println();
        System.out.println("--------------------------");
    }

    /**
     *
     */
    @Test
    public void testConvertUnconvert() {
        System.out.println("Convert & unconvert " + this.source.toString() + " <--> " + this.target.toString() +
                           " on the whole sphere");
        double deltaLongitudeMax = 0.0;
        double deltaLatitudeMax = 0.0;

        for (int latitude = -90; latitude <= 90; latitude++) {
            for (int longitude = 0; longitude < 360; longitude++) {
                final SkyPosition position = source.convertTo(target, longitude, latitude);
                final SkyPosition sourcePosition =
                        target.convertTo(source, position.getLongitude(), position.getLatitude());

                double deltaLongitude = Math.abs(longitude - sourcePosition.getLongitude());
                if (deltaLongitude > 180) {
                    deltaLongitude = 360 - deltaLongitude;
                }

                if (Math.abs(latitude) != 90 && deltaLongitude > deltaLongitudeMax) {
                    deltaLongitudeMax = deltaLongitude;
                }

                final double deltaLatitude = Math.abs(latitude - sourcePosition.getLatitude());
                if (deltaLatitude > deltaLatitudeMax) {
                    deltaLatitudeMax = deltaLatitude;
                }

                if (deltaLatitude > tolerance) {
                    System.out.printf("longitude = %d lat = %d\n", longitude, latitude);
                    System.out.printf("convert: longitude = %20.15f lat = %20.15f\n", position.getLongitude(),
                                      position.getLatitude());
                    System.out.printf("Unconvert : longitude = %20.15f latitude = %20.15f\n",
                                      sourcePosition.getLongitude(), sourcePosition.getLatitude());
                    System.out.println();
                } else if (Math.abs(latitude) != 90 && deltaLongitude > tolerance) {
                    System.out.printf("longitude = %d lat = %d\n", longitude, latitude);
                    System.out.printf("convert: longitude = %20.15f lat = %20.15f\n", position.getLongitude(),
                                      position.getLatitude());
                    System.out.printf("Unconvert : longitude = %20.15f latitude = %20.15f\n",
                                      sourcePosition.getLongitude(), sourcePosition.getLatitude());
                    System.out.println();
                }
            }
        }
        System.out.printf("  Maximum residual (sky): lng: %12.6e  lat: %12.6e\n", deltaLongitudeMax, deltaLatitudeMax);
        assertEquals(1e-12, deltaLongitudeMax, tolerance);
        assertEquals(1e-12, deltaLatitudeMax, tolerance);
    }

}
