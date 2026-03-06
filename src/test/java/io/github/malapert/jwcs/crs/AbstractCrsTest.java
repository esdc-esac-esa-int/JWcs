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
import io.github.malapert.jwcs.datum.FK4;
import io.github.malapert.jwcs.datum.ICRS;
import io.github.malapert.jwcs.datum.FK5;
import io.github.malapert.jwcs.datum.CoordinateReferenceFrame;
import io.github.malapert.jwcs.datum.FK4NoEterms;
import io.github.malapert.jwcs.AbstractJWcs;
import io.github.malapert.jwcs.JWcsFits;
import static io.github.malapert.jwcs.crs.AbstractCrs.convertMatrixEpoch12Epoch2;
import static io.github.malapert.jwcs.crs.AbstractCrs.convertMatrixEqB19502Gal;
import io.github.malapert.jwcs.proj.exception.JWcsError;
import io.github.malapert.jwcs.proj.exception.JWcsException;
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.createRealIdentityMatrix;
import java.io.IOException;
import java.net.URL;
import java.util.logging.Level;
import java.util.logging.Logger;
import nom.tam.fits.Fits;
import nom.tam.fits.FitsException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.apache.commons.math3.linear.RealMatrix;

/**
 *
 * @author Jean-Christophe Malapert
 */
public class AbstractCrsTest {

    private static final double EPSILON_SINGLE = 1e-12;

    /**
     *
     */
    public AbstractCrsTest() {
        //do nothing
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

    @Test
    public void testToString() {
        System.out.println("Test toString");
        AbstractCrs crs1 = new Galactic();
        AbstractCrs crs2 = new Equatorial(new FK5("J2001"));
        String result1 = crs1.toString();
        String result2 = crs2.toString();
        String expectedResult1 = "GALACTIC";
        String expectedResult2 = "EQUATORIAL(FK5(J2001.0))";
        assertEquals(expectedResult2, result2);
        assertEquals(expectedResult1, result1);
        
    }
    
    /**
     *
     */
    @Test
    public void separation() {
        System.out.println("separation");
        final CoordinateReferenceFrame icrs = new ICRS();
        final AbstractCrs sysEqIcrs = new Equatorial(icrs);
        final CoordinateReferenceFrame fk5 = new FK5();
        final AbstractCrs sysEqFK5 = new Equatorial(fk5);
        final SkyPosition pos1 = new SkyPosition(10, 9, sysEqIcrs);
        final SkyPosition pos2 = new SkyPosition(11, 10, sysEqFK5);
        final double separation = SkyPosition.separation(pos1, pos2);
        final double expectedSeparation = 1.4045335865d;
        assertEquals(expectedSeparation, separation, 1e-8);
    }

    /**
     *
     */
    @Test
    public void testConvertFK4B1950ToFK5J2000() {
        System.out.println("convert FK4 B1950 to FK5 J2000");
        final CoordinateReferenceFrame fk4 = new FK4("B1950");
        final AbstractCrs sysEqFK4 = new Equatorial(fk4);
        final CoordinateReferenceFrame fk5 = new FK5("J2000");
        final AbstractCrs sysEqFK5 = new Equatorial(fk5);
        final SkyPosition pos = sysEqFK4.convertTo(sysEqFK5, 0.0d, 0.0d);
        final double expectedLongitude = 0.640691d;
        final double expectedLatitude = 0.27840944d;
        assertEquals(expectedLongitude, pos.getLongitude(), 0.000001);
        assertEquals(expectedLatitude, pos.getLatitude(), 0.000001);
    }

    /**
     *
     */
    @Test
    public void testConvertInverseFK4B1950ToFK5J2000() {
        System.out.println("convert FK4 B1950 to FK5 J2000 and inverse");
        final CoordinateReferenceFrame fk4 = new FK4("B1950");
        final AbstractCrs sysEqFK4 = new Equatorial(fk4);
        final CoordinateReferenceFrame fk5 = new FK5("J2000");
        final AbstractCrs sysEqFK5 = new Equatorial(fk5);
        SkyPosition pos = sysEqFK4.convertTo(sysEqFK5, 30.031d, 10.031d);
        pos = sysEqFK5.convertTo(sysEqFK4, pos.getLongitude(), pos.getLatitude());
        final double expectedLongitude = 30.031d;
        final double expectedLatitude = 10.031d;
        assertEquals(expectedLongitude, pos.getLongitude(), 0.0001);
        assertEquals(expectedLatitude, pos.getLatitude(), 0.0001);
    }

    /**
     *
     */
    @Test
    /**
     * Test based on http://docs.astropy.org/en/stable/coordinates/
     */
    public void testConvertIcrsToGal() {
        System.out.println("convert ICRS to Gal");
        final CoordinateReferenceFrame icrs = new ICRS();
        final AbstractCrs sysEqIcrs = new Equatorial(icrs);
        final AbstractCrs galactic = new Galactic();
        final SkyPosition posInGal = sysEqIcrs.convertTo(galactic, 10.68458d, 41.26917d);
        final double expectedLongitude = 121.174241811d;
        final double expectedLatitude = -21.5728855724d;
        assertEquals(expectedLongitude, posInGal.getLongitude(), 1e-8);
        assertEquals(expectedLatitude, posInGal.getLatitude(), 1e-8);
    }

    /**
     *
     */
    @Test
    /**
     * Test based on http://docs.astropy.org/en/stable/coordinates/
     */
    public void testConvertInverseIcrsToGal() {
        System.out.println("convert ICRS to Gal and Inverse");
        final CoordinateReferenceFrame icrs = new ICRS();
        final AbstractCrs sysEqIcrs = new Equatorial(icrs);
        final AbstractCrs galactic = new Galactic();
        SkyPosition posInGal = sysEqIcrs.convertTo(galactic, 10.68458d, 41.26917d);
        posInGal = galactic.convertTo(sysEqIcrs, posInGal.getLongitude(), posInGal.getLatitude());
        final double expectedLongitude = 10.68458d;
        final double expectedLatitude = 41.26917d;
        assertEquals(expectedLongitude, posInGal.getLongitude(), EPSILON_SINGLE);
        assertEquals(expectedLatitude, posInGal.getLatitude(), EPSILON_SINGLE);
    }

    /**
     *
     */
    @Test
    /**
     * Test based on http://docs.astropy.org/en/stable/coordinates/
     */
    public void testConvertArrayIcrsToGal() {
        System.out.println("convert Array ICRS to Gal");
        final CoordinateReferenceFrame icrs = new ICRS();
        final AbstractCrs sysEqIcrs = new Equatorial(icrs);
        final AbstractCrs galactic = new Galactic();
        SkyPosition posInGal = sysEqIcrs.convertTo(galactic, 10.68458d, 41.26917d);
        final double expectedLongitude1 = posInGal.getLongitude();
        final double expectedLatitude1 = posInGal.getLatitude();

        posInGal = sysEqIcrs.convertTo(galactic, 0.68458d, 1.26917d);
        final double expectedLongitude2 = posInGal.getLongitude();
        final double expectedLatitude2 = posInGal.getLatitude();

        final SkyPosition[] posInGalArray = sysEqIcrs.convertTo(galactic, new double[]{10.68458d, 41.26917d, 0.68458d, 1.26917d});
        assertEquals(expectedLongitude1, posInGalArray[0].getLongitude(), EPSILON_SINGLE);
        assertEquals(expectedLatitude1, posInGalArray[0].getLatitude(), EPSILON_SINGLE);
        assertEquals(expectedLongitude2, posInGalArray[1].getLongitude(), EPSILON_SINGLE);
        assertEquals(expectedLatitude2, posInGalArray[1].getLatitude(), EPSILON_SINGLE);
    }

    /**
     *
     */
    @Test
    /**
     * Test based on http://docs.astropy.org/en/stable/coordinates/
     */
    public void testConvertIcrsToFK5J2000() {
        System.out.println("convert ICRS to FK5 J2000");
        final CoordinateReferenceFrame icrs = new ICRS();
        final AbstractCrs sysEqIcrs = new Equatorial(icrs);
        final CoordinateReferenceFrame fk5 = new FK5("J2000");
        final AbstractCrs eqFk5 = new Equatorial(fk5);
        final SkyPosition posInFk5 = sysEqIcrs.convertTo(eqFk5, 10.68458d, 41.26917d);
        final double expectedLongitude = 10.6845915393d;
        final double expectedLatitude = 41.2691714591d;
        assertEquals(expectedLongitude, posInFk5.getLongitude(), 1e-10);
        assertEquals(expectedLatitude, posInFk5.getLatitude(), 1e-10);
    }

    /**
     *
     */
    @Test
    /**
     * Test based on http://docs.astropy.org/en/stable/coordinates/
     */
    public void testConvertInverseIcrsToFK5J2000() {
        System.out.println("convert ICRS to FK5 J2000 and inverse");
        final CoordinateReferenceFrame icrs = new ICRS();
        final AbstractCrs sysEqIcrs = new Equatorial(icrs);
        final CoordinateReferenceFrame fk5 = new FK5("J2000");
        final AbstractCrs eqFk5 = new Equatorial(fk5);
        SkyPosition posInFk5 = sysEqIcrs.convertTo(eqFk5, 10.68458d, 41.26917d);
        posInFk5 = eqFk5.convertTo(sysEqIcrs, posInFk5.getLongitude(), posInFk5.getLatitude());
        final double expectedLongitude = 10.68458d;
        final double expectedLatitude = 41.26917d;
        assertEquals(expectedLongitude, posInFk5.getLongitude(), EPSILON_SINGLE);
        assertEquals(expectedLatitude, posInFk5.getLatitude(), EPSILON_SINGLE);
    }

    /**
     *
     */
    @Test
    /**
     * Test based on http://docs.astropy.org/en/stable/coordinates/
     */
    public void testConvertToFK5J1975() {
        System.out.println("convert FK5J2000 to FK5 J1975");
        final CoordinateReferenceFrame icrs = new ICRS();
        final AbstractCrs sysEqIcrs = new Equatorial(icrs);
        final CoordinateReferenceFrame fk5 = new FK5("J2000");
        final AbstractCrs eqFk5J2000 = new Equatorial(fk5);
        final SkyPosition posInFk5J2000 = sysEqIcrs.convertTo(eqFk5J2000, 10.68458d, 41.26917d);
        final CoordinateReferenceFrame fk5_1975 = new FK5("J1975");
        final AbstractCrs eqFk5J1975 = new Equatorial(fk5_1975);
        final SkyPosition posInFk5J1975 = eqFk5J2000.convertTo(eqFk5J1975, posInFk5J2000.getLongitude(), posInFk5J2000.getLatitude());
        final double expectedLongitude = 10.3420913461d;
        final double expectedLatitude = 41.1323211229d;
        assertEquals(expectedLongitude, posInFk5J1975.getLongitude(), 0.00002);
        assertEquals(expectedLatitude, posInFk5J1975.getLatitude(), 0.00002);
    }

    /**
     *
     */
    @Test
    /**
     * Test based on http://docs.astropy.org/en/stable/coordinates/
     */
    public void testConvertInverseToFK5J1975() {
        System.out.println("convert FK5J2000 to FK5 J1975 and inverse");

        final CoordinateReferenceFrame fk5 = new FK5("J2000");
        final AbstractCrs eqFk5J2000 = new Equatorial(fk5);

        final CoordinateReferenceFrame fk51975 = new FK5("J1975");
        final AbstractCrs eqFk5J1975 = new Equatorial(fk51975);

        SkyPosition posInFk5J1975 = eqFk5J2000.convertTo(eqFk5J1975, 10.3420913461d, 41.1323211229d);
        posInFk5J1975 = eqFk5J1975.convertTo(eqFk5J2000, posInFk5J1975.getLongitude(), posInFk5J1975.getLatitude());
        final double expectedLongitude = 10.3420913461d;
        final double expectedLatitude = 41.1323211229d;
        assertEquals(expectedLongitude, posInFk5J1975.getLongitude(), 1e-8);
        assertEquals(expectedLatitude, posInFk5J1975.getLatitude(), 1e-8);
    }

    /**
     * Test of convertTo method, of class AbstractCrs.
     */
    @Test
    public void testIcrsToFK5() {
        System.out.println("convert ICRS To FK5");
        final CoordinateReferenceFrame ref = new ICRS();
        final AbstractCrs sys1 = new Equatorial(ref);
        final AbstractCrs sys2 = new Equatorial(new FK5());
        final SkyPosition position = sys1.convertTo(sys2, 182.63867d, 39.401167d);
        final double expectedLongitude = 182.63867d;
        final double expectedLatitude = 39.401165d;
        assertEquals(expectedLongitude, position.getLongitude(), 0.00001);
        assertEquals(expectedLatitude, position.getLatitude(), 0.000001);
    }

    /**
     *
     */
    @Test
    public void testIcrsToFK5Inverse() {
        System.out.println("convert ICRS To FK5 and Inverse");
        final CoordinateReferenceFrame ref = new ICRS();
        final AbstractCrs sys1 = new Equatorial(ref);
        final AbstractCrs sys2 = new Equatorial(new FK5());
        SkyPosition position = sys1.convertTo(sys2, 182.63867d, 39.401167d);
        final double expectedLongitude = 182.63867d;
        final double expectedLatitude = 39.401167d;
        position = sys2.convertTo(sys1, position.getLongitude(), position.getLatitude());
        assertEquals(expectedLongitude, position.getLongitude(), EPSILON_SINGLE);
        assertEquals(expectedLatitude, position.getLatitude(), EPSILON_SINGLE);
    }

    /**
     *
     */
    @Test
    public void testIcrsToFK4() {
        System.out.println("convert ICRS To FK4");
        final CoordinateReferenceFrame ref = new ICRS();
        final AbstractCrs sys1 = new Equatorial(ref);
        final AbstractCrs sys2 = new Equatorial(new FK4());
        final SkyPosition position = sys1.convertTo(sys2, 182.63867, 39.401167);
        final double expectedLongitude = 182.0073;
        final double expectedLatitude = 39.679217;
        assertEquals(expectedLongitude, position.getLongitude(), 0.0001);
        assertEquals(expectedLatitude, position.getLatitude(), 0.0001);
    }

    /**
     *
     */
    @Test
    public void testIcrsToFK4Inverse() {
        System.out.println("convert ICRS To FK4 and inverse");
        final CoordinateReferenceFrame ref = new ICRS();
        final AbstractCrs sys1 = new Equatorial(ref);
        final AbstractCrs sys2 = new Equatorial(new FK4());
        SkyPosition position = sys1.convertTo(sys2, 182.63867d, 39.401167d);
        position = sys2.convertTo(sys1, position.getLongitude(), position.getLatitude());
        final double expectedLongitude = 182.63867d;
        final double expectedLatitude = 39.401167d;
        assertEquals(expectedLongitude, position.getLongitude(), 0.0001);
        assertEquals(expectedLatitude, position.getLatitude(), 0.0001);
    }

    /**
     *
     */
    @Test
    public void testIcrsToGalactic() {
        System.out.println("convert ICRS To Galactic");
        final CoordinateReferenceFrame ref = new ICRS();
        final AbstractCrs sys1 = new Equatorial(ref);
        final AbstractCrs sys2 = new Galactic();
        final SkyPosition position = sys1.convertTo(sys2, 182.63867, 39.401167);
        final double expectedLongitude = 155.08125;
        final double expectedLatitude = 75.068157;
        assertEquals(expectedLongitude, position.getLongitude(), 0.00001);
        assertEquals(expectedLatitude, position.getLatitude(), 0.00001);
    }

    /**
     *
     */
    @Test
    public void testIcrsToGalacticInverse() {
        System.out.println("convert ICRS To Galactic and inverse");
        final CoordinateReferenceFrame ref = new ICRS();
        final AbstractCrs sys1 = new Equatorial(ref);
        final AbstractCrs sys2 = new Galactic();
        SkyPosition position = sys1.convertTo(sys2, 182.63867d, 39.401167d);
        position = sys2.convertTo(sys1, position.getLongitude(), position.getLatitude());
        final double expectedLongitude = 182.63867d;
        final double expectedLatitude = 39.401167d;
        assertEquals(expectedLongitude, position.getLongitude(), EPSILON_SINGLE);
        assertEquals(expectedLatitude, position.getLatitude(), EPSILON_SINGLE);
    }

    /**
     *
     */
    @Test
    public void testConvertDegreeToHexa() {
        System.out.println("convert Degrees to hms/dms");
        final SkyPosition pos = new SkyPosition(182.63867, 39.401167, new Equatorial());
        assertEquals("12:10:33.281", pos.getLongitudeAsSexagesimal());
        assertEquals("+39:24:04.20", pos.getLatitudeAsSexagesimal());
    }

    /**
     *
     */
    @Test
    public void testICRS2GAL() {
        System.out.println("ICRS <--> GAL");
        final double[][] coordinates = new double[][]{
            {0, 0},
            {180, 60},
            {359, 60},
            {86, -35}
        };
        final AbstractCrs sk1 = new Equatorial();
        final AbstractCrs sk2 = new Galactic();
        for (final double[] coordinate : coordinates) {
            SkyPosition[] result = sk1.convertTo(sk2, coordinate);
            result = sk2.convertTo(sk1, result[0].getDoubleArray());
            assertArrayEquals(coordinate, result[0].getDoubleArray(), EPSILON_SINGLE);
        }
    }

    /**
     *
     */
    @Test
    public void testICRS2SUPERGAL() {
        System.out.println("ICRS <--> SUPERGAL");
        final double[][] coordinates = new double[][]{
            {0, 0},
            {180, 60},
            {359, 60},
            {86, -35}
        };
        final AbstractCrs sk1 = new Equatorial();
        final AbstractCrs sk2 = new SuperGalactic();
        for (final double[] coordinate : coordinates) {
            SkyPosition[] result = sk1.convertTo(sk2, coordinate);
            result = sk2.convertTo(sk1, result[0].getDoubleArray());
            assertArrayEquals(coordinate, result[0].getDoubleArray(), EPSILON_SINGLE);
        }
    }

    /**
     *
     */
    @Test
    public void testGAL2SUPERGAL() {
        System.out.println("GAL <--> SUPERGAL");
        final double[][] coordinates = new double[][]{
            {0, 0},
            {180, 60},
            {359, 60},
            {86, -35}
        };
        final AbstractCrs sk1 = new Galactic();
        final AbstractCrs sk2 = new SuperGalactic();
        for (final double[] coordinate : coordinates) {
            SkyPosition[] result = sk1.convertTo(sk2, coordinate);
            result = sk2.convertTo(sk1, result[0].getDoubleArray());
            assertArrayEquals(coordinate, result[0].getDoubleArray(), EPSILON_SINGLE);
        }
    }

    /**
     *
     */
    @Test
    public void testICRS2FK5() {
        System.out.println("ICRS <--> FK5");
        final double[][] coordinates = new double[][]{
            {0, 0},
            {180, 60},
            {359, 60},
            {86, -35}
        };
        final AbstractCrs sk1 = new Equatorial();
        final AbstractCrs sk2 = new Equatorial(new FK5());
        for (final double[] coordinate : coordinates) {
            SkyPosition[] result = sk1.convertTo(sk2, coordinate);
            result = sk2.convertTo(sk1, result[0].getDoubleArray());
            assertArrayEquals(coordinate, result[0].getDoubleArray(), EPSILON_SINGLE);
        }
    }

    /**
     *
     */
    @Test
    public void testFK5J20002FK5J1950() {
        System.out.println("FK5(J2000) <--> FK5(1950)");
        final double[][] coordinates = new double[][]{
            {1, 0},
            {180, 60},
            {359, 60},
            {86, -35}
        };
        final AbstractCrs sk1 = new Equatorial(new FK5("J2000"));
        final AbstractCrs sk2 = new Equatorial(new FK5("B1950"));
        for (final double[] coordinate : coordinates) {
            SkyPosition[] result = sk1.convertTo(sk2, coordinate);
            result = sk2.convertTo(sk1, result[0].getDoubleArray());
            assertArrayEquals(coordinate, result[0].getDoubleArray(), 1e-7);
        }
    }

    /**
     *
     */
    @Test
    public void testFK5J20002FK4() {
        System.out.println("FK5(J2000) <--> FK4");
        final double[][] coordinates = new double[][]{
            {1, 0},
            {180, 60},
            {359, 60},
            {86, -35}
        };
        final AbstractCrs sk1 = new Equatorial(new FK5("J2000"));
        final AbstractCrs sk2 = new Equatorial(new FK4());
        for (final double[] coordinate : coordinates) {
            SkyPosition[] result = sk1.convertTo(sk2, coordinate);
            result = sk2.convertTo(sk1, result[0].getDoubleArray());
            assertArrayEquals(coordinate, result[0].getDoubleArray(), 1e-10);
        }
    }

    /**
     *
     */
    @Test
    public void testCurrentSkySystemToGalacticFromFITS() {
        System.out.println("From Fits <--> Gal");
        try {

            URL url = getClass().getClassLoader().getResource("WFPC2ASSNu5780205bx.fits");
            assertNotNull(url);
            final AbstractJWcs wcs = new JWcsFits(new Fits(url));
            wcs.doInit();
            // convert pixel(1,1) to Sky
            final double[] posOrigin = wcs.pix2wcs(1, 1);
            //convert (ra,dec) To galactic
            final AbstractCrs sysOrigin = wcs.getCrs();
            final AbstractCrs sysTarget = new Galactic();
            assertEquals("EQUATORIAL", sysOrigin.getCoordinateReferenceSystem().name());
            final SkyPosition skyPosTarget = sysOrigin.convertTo(sysTarget, posOrigin[0], posOrigin[1]);
            
            //convert skyPos (galactic frame) to Equatorial
            final SkyPosition newPosOrigin = sysTarget.convertTo(sysOrigin, skyPosTarget.getLongitude(), skyPosTarget.getLatitude());
            assertArrayEquals(posOrigin, newPosOrigin.getDoubleArray(), EPSILON_SINGLE);
            // convert newPos to camera
            final double[] returnedPosCamera = wcs.wcs2pix(newPosOrigin.getLongitude(), newPosOrigin.getLatitude());
            assertArrayEquals(new double[]{1, 1}, returnedPosCamera, 0.5);//Set the precision at the half pixel
        } catch (FitsException | IOException | JWcsException ex) {
            Logger.getLogger(AbstractCrsTest.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
    
    /**
     *
     */
    @Test
    public void testSkyMatrixFK4toFK5() {
        final CoordinateReferenceFrame fk4 = new FK4();
        final CoordinateReferenceFrame fk5 = new FK5();
        final AbstractCrs eq1 = new Equatorial(fk4);
        final AbstractCrs eq2 = new Equatorial(fk5);
        final RealMatrix rotationMatrix = eq1.getRotationMatrix(eq2);
        final RealMatrix etermsIn = AbstractCrs.getEterms(eq1);
        final RealMatrix etermsOut = AbstractCrs.getEterms(eq2);
        assertArrayEquals(new double[]{9.99925679e-01,  -1.11814832e-02,  -4.85900382e-03}, rotationMatrix.getRow(0), 1e-9);
        assertArrayEquals(new double[]{1.11814832e-02,   9.99937485e-01,  -2.71625947e-05}, rotationMatrix.getRow(1), 1e-9);
        assertArrayEquals(new double[]{4.85900377e-03,  -2.71702937e-05,   9.99988195e-01}, rotationMatrix.getRow(2), 1e-9);
        assertArrayEquals(new double[]{-1.6255503575995309e-06,-3.1918587795578522e-07,-1.3842701121066153e-07}, etermsIn.getRow(0),1e-20);
        assertNull(etermsOut);
    }
    
    /**
     *
     */
    @Test
    public void testSkyMatrixFK4NOEB1950toFK4B1950() {
        final CoordinateReferenceFrame fk4 = new FK4("B1950");
        final CoordinateReferenceFrame fk4NOE = new FK4NoEterms("B1950");
        final AbstractCrs eq1 = new Equatorial(fk4NOE);
        final AbstractCrs eq2 = new Equatorial(fk4);
        final RealMatrix rotationMatrix = eq1.getRotationMatrix(eq2);
        final RealMatrix etermsIn = AbstractCrs.getEterms(eq1);
        final RealMatrix etermsOut = AbstractCrs.getEterms(eq2);
        assertArrayEquals(new double[]{1.,  0.,  -0.}, rotationMatrix.getRow(0), 1e-9);
        assertArrayEquals(new double[]{0.,   1.,  0.}, rotationMatrix.getRow(1), 1e-9);
        assertArrayEquals(new double[]{0.,  0.,   1.}, rotationMatrix.getRow(2), 1e-9);
        assertArrayEquals(new double[]{-1.6255503575995309e-06,-3.1918587795578522e-07,-1.3842701121066153e-07}, etermsOut.getRow(0),1e-20);
        assertNull(etermsIn);
    }    
    
    /**
     *
     */
    @Test
    public void testSkyMatrixFK4B1950J19835toFK52000() {
        final CoordinateReferenceFrame fk4 = new FK4("B1950", "J1983.5");
        final CoordinateReferenceFrame fk5 = new FK5("J2000");
        final AbstractCrs eq1 = new Equatorial(fk4);
        final AbstractCrs eq2 = new Equatorial(fk5);
        final RealMatrix rotationMatrix = eq1.getRotationMatrix(eq2);
        final RealMatrix etermsIn = AbstractCrs.getEterms(eq1);
        final RealMatrix etermsOut = AbstractCrs.getEterms(eq2);
        assertArrayEquals(new double[]{9.99925679e-01,  -1.11818698e-02,  -4.85829658e-03}, rotationMatrix.getRow(0), 1e-9);
        assertArrayEquals(new double[]{1.11818699e-02,   9.99937481e-01,  -2.71546879e-05}, rotationMatrix.getRow(1), 1e-9);
        assertArrayEquals(new double[]{4.85829648e-03,  -2.71721706e-05,   9.99988198e-01}, rotationMatrix.getRow(2), 1e-9);
        assertArrayEquals(new double[]{-1.6255503575995309e-06,-3.1918587795578522e-07,-1.3842701121066153e-07}, etermsIn.getRow(0),1e-20);
        assertNull(etermsOut);
    } 
    
    /**
     * Test of convertMatrixEqB19502Gal method, of class Utility.
     */
    @Test
    public void testMatrixEqB19502Gal() {
        System.out.println("MatrixEqB19502Gal");
        final RealMatrix result1 = convertMatrixEqB19502Gal();
        final RealMatrix result2 = convertMatrixEqB19502Gal().transpose();
        final RealMatrix result = result1.multiply(result2);
        final RealMatrix expResult = createRealIdentityMatrix(3);
        assertArrayEquals(expResult.getRow(0), result.getRow(0), NumericalUtility.DOUBLE_TOLERANCE);
        assertArrayEquals(expResult.getRow(1), result.getRow(1), NumericalUtility.DOUBLE_TOLERANCE);
        assertArrayEquals(expResult.getRow(2), result.getRow(2), NumericalUtility.DOUBLE_TOLERANCE);
    }

    /**
     * Test of convertMatrixEpoch12Epoch2 method, of class Utility.
     */
    @Test
    public void testMatrixEpoch12Epoch2() {
        System.out.println("MatrixEpoch12Epoch2");
        final double epoch1 = 1940d;
        final double epoch2 = 1960d;
        final CoordinateReferenceFrame.ReferenceFrame s1 = CoordinateReferenceFrame.ReferenceFrame.FK4;
        final CoordinateReferenceFrame.ReferenceFrame s2 = CoordinateReferenceFrame.ReferenceFrame.FK5;
        final double epobs = 1950d;
        final RealMatrix result = convertMatrixEpoch12Epoch2(epoch1, epoch2, s1, s2, epobs);
        assertArrayEquals(new double[]{9.99988107e-01,  -4.47301372e-03,  -1.94362889e-03}, result.getRow(0), 1e-9);
        assertArrayEquals(new double[]{4.47301372e-03,   9.99989996e-01,  -4.34712255e-06}, result.getRow(1), 1e-9);
        assertArrayEquals(new double[]{1.94362889e-03,  -4.34680782e-06,   9.99998111e-01}, result.getRow(2), 1e-9);
    }    
    
    @Test
    public void testExceptionConvertTo() {
        System.out.println("Exception in convertTo");
        double[] longlat = new double[]{0,0,1,6,2};
        AbstractCrs crs = new Galactic();
        final JWcsError expected = new JWcsError("coordinates should be an array containing a set of [longitude, latitude]");
        JWcsError result = null;
        try {
            crs.convertTo(new SuperGalactic(), longlat);
        } catch (JWcsError error) {
            result = error;
        } finally {
            assertEquals(expected.toString(), result.toString());
        }
    }
    
    @Test
    public void testExceptionCoordinates() {
        System.out.println("Exception in convertTo");
        AbstractCrs crs = new Galactic();        
        JWcsError expected = new JWcsError("longitude must be in [0,360]");
        JWcsError result = null;
        try {
            crs.convertTo(new SuperGalactic(), 400, -90);
        } catch (JWcsError error) {
            result = error;
        } finally {
            assertEquals(expected.toString(), result.toString());
        }
        
        expected = new JWcsError("latitude must be in [-90,90]");
        result = null;
        try {
            crs.convertTo(new SuperGalactic(), 360, -91);
        } catch (JWcsError error) {
            result = error;
        } finally {
            assertEquals(expected.toString(), result.toString());
        }  
        
        expected = new JWcsError("longitude must be in [0,360] and latitude in [-90,90]");
        result = null;
        try {
            crs.convertTo(new SuperGalactic(), 400, -91);
        } catch (JWcsError error) {
            result = error;
        } finally {
            assertEquals(expected.toString(), result.toString());
        }         
    }    
    
}
