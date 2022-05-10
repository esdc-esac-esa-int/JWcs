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

import io.github.malapert.jwcs.crs.Ecliptic;
import io.github.malapert.jwcs.crs.Equatorial;
import io.github.malapert.jwcs.datum.FK4;
import io.github.malapert.jwcs.datum.FK4NoEterms;
import io.github.malapert.jwcs.datum.FK5;
import io.github.malapert.jwcs.crs.Galactic;
import io.github.malapert.jwcs.datum.ICRS;
import io.github.malapert.jwcs.crs.AbstractCrs;
import io.github.malapert.jwcs.proj.AbstractProjection;
import io.github.malapert.jwcs.proj.AbstractProjection.ProjectionParameter;
import io.github.malapert.jwcs.proj.SZP;
import io.github.malapert.jwcs.proj.TPV;
import io.github.malapert.jwcs.proj.ZPN;
import io.github.malapert.jwcs.proj.exception.BadProjectionParameterException;
import io.github.malapert.jwcs.proj.exception.JWcsException;
import io.github.malapert.jwcs.proj.exception.JWcsError;
import io.github.malapert.jwcs.proj.exception.ProjectionException;
import io.github.malapert.jwcs.utility.TimeUtility;
import java.text.ParseException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import nom.tam.fits.Header;
import nom.tam.fits.HeaderCard;
import io.github.malapert.jwcs.datum.CoordinateReferenceFrame;
import io.github.malapert.jwcs.proj.BON;
import io.github.malapert.jwcs.proj.CEA;
import io.github.malapert.jwcs.proj.COD;
import io.github.malapert.jwcs.proj.COE;
import io.github.malapert.jwcs.proj.COO;
import io.github.malapert.jwcs.proj.COP;
import io.github.malapert.jwcs.proj.SIN;
import static io.github.malapert.jwcs.utility.NumericalUtility.createRealMatrix;
import static io.github.malapert.jwcs.utility.NumericalUtility.inverse;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Locale;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

/**
 * The FITS "World Coordinate System" (WCS) standard defines keywords and usage
 * that provide for the description of astronomical coordinate systems in a FITS
 * image header.
 *
 * <p>
 * This class provides methods to compute WCS transformation.
 * </p>
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public abstract class AbstractJWcs implements JWcsKeyProvider {

    /**
     * Maximum longitude value in degrees.
     */
    public static final int MAX_LONGITUDE = 360;

    /**
     * Minimum longitude value in degrees.
     */
    public static final int MIN_LONGITUDE = 0;

    /**
     * Minimum latitude value in degrees.
     */
    public static final int MIN_LATITUDE = -90;

    /**
     * Maximum latitude value in degrees.
     */
    public static final int MAX_LATITUDE = 90;

    /**
     * Number of axes. 
     * 
     * <p>2 for an image
     */
    public static final String NAXIS = "NAXIS";
    /**
     * Number of pixels along X axis.
     */
    public static final String NAXIS1 = "NAXIS1";
    /**
     * Number of pixels along Y axis.
     */
    public static final String NAXIS2 = "NAXIS2";
    /**
     * Reference along X axis in pixel frame. 
     * 
     * <p>This keyword is required for projection computation.
     */
    public static final String CRPIX1 = "CRPIX1";
    /**
     * Reference along Y axis in pixel frame. 
     * 
     * <p>This keyword is required for projection computation.
     */
    public static final String CRPIX2 = "CRPIX2";
    /**
     * Reference along longitude in degrees in celestial frame. 
     * 
     * <p>This keyword is required for projection computation.
     */
    public static final String CRVAL1 = "CRVAL1";
    /**
     * Reference along latitude in degrees in celestial frame. 
     * 
     * <p>This keyword is required for projection computation.
     */
    public static final String CRVAL2 = "CRVAL2";
    /**
     * AbstractProjection type along X axis. 
     * 
     * <p>This keyword is required for projection computation.
     */
    public static final String CTYPE1 = "CTYPE1";
    /**
     * AbstractProjection type along Y axis. 
     * 
     * <p>This keyword is required for projection computation.
     */
    public static final String CTYPE2 = "CTYPE2";
    /**
     * Scale (degrees / pixel) and rotation matrix. 
     * 
     * <p>For projection computation, information about scale and rotation are 
     * needed. Either the CD matrix is provided or the following element 
     * (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1, CDELT2).
     */
    public static final String CD11 = "CD1_1";
    /**
     * Scale (degrees / pixel) and rotation matrix. 
     * 
     * <p>For projection computation, information about scale and rotation are
     * needed. Either the CD matrix is provided or the following element 
     * (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1, CDELT2).
     */
    public static final String CD12 = "CD1_2";
    /**
     * Scale (degrees / pixel) and rotation matrix. 
     * 
     * <p>For projection computation, information about scale and rotation are 
     * needed. Either the CD matrix is provided or the following element 
     * (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1, CDELT2).
     */
    public static final String CD21 = "CD2_1";
    /**
     * Scale (degrees / pixel) and rotation matrix. 
     * 
     * <p>For projection computation, information about scale and rotation are
     * needed. Either the CD matrix is provided or the following element 
     * (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1, CDELT2).
     */
    public static final String CD22 = "CD2_2";
    /**
     * Unit along X axis.
     */
    public static final String CUNIT1 = "CUNIT1";
    /**
     * Unit along Y axis.
     */
    public static final String CUNIT2 = "CUNIT2";
    /**
     * Scale (degrees / pixel) along X axis when CD matrix is not defined. 
     * 
     * <p>For projection computation, information about scale and rotation are needed.
     * Either the CD matrix is provided or the following element (CDELT1,
     * CDELT2, CROTA2) or (PC matrix, CDELT1, CDELT2).
     */
    public static final String CDELT1 = "CDELT1";
    /**
     * Scale (degrees / pixel) along X axis when CD matrix is not defined. 
     * 
     * <p>For projection computation, information about scale and rotation are needed.
     * Either the CD matrix is provided or the following element (CDELT1,
     * CDELT2, CROTA2) or (PC matrix, CDELT1, CDELT2).
     */
    public static final String CDELT2 = "CDELT2";
    /**
     * For projection computation, information about scale and rotation are
     * needed. 
     * 
     * <p>Either the CD matrix is provided or the following element
     * (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1, CDELT2).
     */
    public static final String CROTA2 = "CROTA2";
    /**
     * Equinox value.
     */
    public static final String EQUINOX = "EQUINOX";
    /**
     * Deformation matrix. 
     * 
     * <p>For projection computation, information about scale
     * and rotation are needed. Either the CD matrix is provided or the
     * following element (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1,
     * CDELT2).
     */
    public static final String PC11 = "PC1_1";
    /**
     * Deformation matrix. 
     * 
     * <p>For projection computation, information about scale
     * and rotation are needed. Either the CD matrix is provided or the
     * following element (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1,
     * CDELT2).
     */
    public static final String PC12 = "PC1_2";
    /**
     * Deformation matrix. 
     * 
     * <p>For projection computation, information about scale
     * and rotation are needed. Either the CD matrix is provided or the
     * following element (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1,
     * CDELT2).
     */
    public static final String PC21 = "PC2_1";
    /**
     * Deformation matrix. 
     * 
     * <p>For projection computation, information about scale
     * and rotation are needed. Either the CD matrix is provided or the
     * following element (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1,
     * CDELT2).
     */
    public static final String PC22 = "PC2_2";
    /**
     * Deformation matrix.
     */
    public static final String PV11 = "PV1_1";
    /**
     * Deformation matrix.
     */
    public static final String PV12 = "PV1_2";
    /**
     * Deformation matrix.
     */
    public static final String PV13 = "PV1_3";
    /**
     * Deformation matrix.
     */
    public static final String PV14 = "PV1_4";
    /**
     * Deformation matrix.
     */
    public static final String PV20 = "PV2_0";
    /**
     * Deformation matrix.
     */
    public static final String PV21 = "PV2_1";
    /**
     * Deformation matrix.
     */
    public static final String PV22 = "PV2_2";
    /**
     * Deformation matrix.
     */
    public static final String PV23 = "PV2_3";
    /**
     * lontpole.
     */
    public static final String LONPOLE = "LONPOLE";
    /**
     * latpole.
     */
    public static final String LATPOLE = "LATPOLE";
    /**
     * Reference system.
     */
    public static final String RADESYS = "RADESYS";

    /**
     * AbstractProjection object.
     */
    private AbstractProjection proj;

    /**
     * CD matrix : scale and rotation.
     */
    private RealMatrix cd;

    /**
     * Inverse CD matrix.
     */
    private RealMatrix cdInverse;

    /**
     * LOG.
     */
    protected static final Logger LOG = Logger.getLogger(AbstractJWcs.class.getName());

    /**
     * Initialize the WCS Object.
     *
     * <p>The WCS object is initialized by doing the following steps:
     * <ul>
     * <li>creates the projection</li>
     * <li>creates the CD matrix</li>
     * <li>creates the CD matrix inverse</li>
     * <li>checks the WCS</li>
     * </ul>
     *
     * @throws JWcsException When WCS is not valid
     */
    protected final void init() throws JWcsException {
        checkWcs();
        setProj(createProjection());
        setCd(createCdMatrix());
        setCdInverse(inverse(getCd()).transpose());
    }

    /**
     * Checks WCS keywords.
     *
     * @throws JWcsException When WCS is not valid
     */
    protected abstract void checkWcs() throws JWcsException;

    /**
     * Checks if the WCS header is valid.
     *
     * @param hdr Header FITS
     * @return True when the Header is valid otherwise False.
     */
    public static boolean isValidWcs(final Header hdr) {
        final AbstractJWcs wcs = new JWcsFits(hdr);
        boolean result;
        try {
            wcs.checkWcs();
            result = true;
        } catch (JWcsException ex) {
            result = false;
        }
        return result;
    }

    /**
     * Returns the modified Julian date.
     * 
     * <p>Returns the value of the MJD-OBS keyword when it is present otherwise
     * returns the value of the DATE-OBS and convert it on the modified Julian
     * date.
     *
     * @return the Modified Julia Date.
     * @throws JWcsError Cannot find or compute MJD-OBS
     */
    private String getMJDObs() {
        String mjd;
        if (hasKeyword("MJD-OBS")) {
            mjd = String.valueOf(getValueAsFloat("MJD-OBS"));
        } else if (hasKeyword("DATE-OBS")) {
            try {
                mjd = String.valueOf(TimeUtility.convertISOToModifiedJulianDate(getValueAsString("DATE-OBS")));
            } catch (ParseException ex) {
                mjd = null;
            }
        } else {
            throw new JWcsError("Cannot find or compute MJD-OBS");
        }
        return mjd;
    }

    /**
     * Returns the reference system.
     * 
     * <p>To find the reference system, the algorithm proceed as it:
     * <ul>
     * <li>Gets the RADESYS value when the keyword is found
     * <li>Otherwise gets the EQUINOX value when the keyword is found and select
     * either FK4 or FK5 according to the equinox value
     * <li>Otherwise ICRS is set
     * </ul>
     *
     * @return the coordinate reference frame
     * @throws JWcsError the coordinate reference frame is not supported
     */
    private CoordinateReferenceFrame getReferenceFrame() {
        CoordinateReferenceFrame refSystem;
        if (hasKeyword(RADESYS)) {
            final String radecsys = getValueAsString(RADESYS);
            refSystem = coordinateReferenceFrameFactory(radecsys);
        } else if (hasKeyword(EQUINOX)) {
            final float equinox = getValueAsFloat(EQUINOX);
            refSystem = coordinateFK4orFK5(equinox);
        } else {
            // RADESYSa defaults to ICRS if both RADESYS and EQUINOX are absents.
            refSystem = new ICRS();
        }
        return refSystem;
    }

    /**
     * Creates a FK4 or FK5 coordinate reference frame depending on the equinox.
     *
     * <p>Create a FK4 coordinate reference frame with a Besselian epoch of equinox
     * and a Modified Julian date as epoch of observation when equinox is
     * smaller then 1984. Otherwise, a FK5 coordinate refrence frame is created
     * based on a Julian date as epoch of observation.
     *
     * @param equinox epoch of the equinox
     * @return a FK4 or FK5 coordinate reference
     */
    private CoordinateReferenceFrame coordinateFK4orFK5(final float equinox) {
        final CoordinateReferenceFrame refSystem;
        if (equinox < 1984.0) {
            refSystem = new FK4("B" + equinox);
            final String mjdObs = getMJDObs();
            if (mjdObs != null && !mjdObs.isEmpty()) {
                refSystem.setEpochObs("MJD" + Float.valueOf(mjdObs));
            }
        } else {
            refSystem = new FK5("J" + equinox);
        }
        return refSystem;
    }

    /**
     * Creates a factory to create the right coordinate reference frame based on
     * the RADESYS keyword value.
     *
     * @param radecSys value of the RADESYS keyword
     * @return a coordinate reference frame
     */
    private CoordinateReferenceFrame coordinateReferenceFrameFactory(final String radecSys) {
        final CoordinateReferenceFrame refSystem;
        switch (radecSys) {
            case "ICRS":
                refSystem = new ICRS();
                break;
            case "FK5":
                refSystem = createFK5();
                break;
            case "FK4":
                refSystem = createFK4();
                break;
            case "FK4-NO-E":
                refSystem = createFK4NOE();
                break;
            default:
                throw new JWcsError("The coordinate reference frame, " + radecSys + " is not supported");
        }
        return refSystem;
    }

    /**
     * Creates a FK5 coordinate reference frame.
     *
     * @return a FK5 coordinate reference frame
     */
    private CoordinateReferenceFrame createFK5() {
        final CoordinateReferenceFrame refSystem = new FK5();
        if (hasKeyword(EQUINOX)) {
            refSystem.setEquinox("J" + getValueAsFloat(EQUINOX));
        }
        return refSystem;
    }

    /**
     * Creates a FK4 coordinate reference frame.
     *
     * @return a FK4 coordinate reference frame
     */
    private CoordinateReferenceFrame createFK4() {
        final CoordinateReferenceFrame refSystem = new FK4();
        if (hasKeyword(EQUINOX)) {
            refSystem.setEquinox("B" + getValueAsFloat(EQUINOX));
        }
        final String mjdObs = getMJDObs();
        if (mjdObs != null && !mjdObs.isEmpty()) {
            refSystem.setEpochObs("MJD" + Float.valueOf(mjdObs));
        }
        return refSystem;
    }

    /**
     * Creates a FK4 coordinate reference frame without Eterms.
     *
     * @return a FK4 coordinate reference frame without Eterms
     */
    private CoordinateReferenceFrame createFK4NOE() {
        final CoordinateReferenceFrame refSystem = new FK4NoEterms();
        if (hasKeyword(EQUINOX)) {
            refSystem.setEquinox("B" + getValueAsFloat(EQUINOX));
        }
        final String mjdObs = getMJDObs();
        if (mjdObs != null && !mjdObs.isEmpty()) {
            refSystem.setEpochObs("MJD" + Float.valueOf(mjdObs));
        }
        return refSystem;
    }

    /**
     * Returns the coordinate reference system.
     * 
     * <p>The coordinate reference system is found according to the CTYPE1 keyword.
     *
     * @return the coordinate reference system
     * @throws JWcsError The coordinate reference system is not supported
     */
    public AbstractCrs getCrs() {
        final AbstractCrs crs;
        if (hasKeyword("CTYPE1")) {
            String ctype1 = getValueAsString("CTYPE1");
            ctype1 = ctype1.substring(0, ctype1.indexOf('-'));
            final CoordinateReferenceFrame refSystem = getReferenceFrame();
            crs = coordinateReferenceSystemFactory(ctype1, refSystem);
        } else {
            throw new JWcsError("cannot find the coordinate reference system.");
        }

        return crs;

    }
    
    /**
     * Create a coordinate reference system based on the coordinate system and
     * the coordinate reference frame.
     *
     * @param ctype keyword value of CTYPE1
     * @param refSystem coordinate reference frame
     * @return a coordinate reference system
     * @throws JWcsError the CRS is not supported
     */
    private AbstractCrs coordinateReferenceSystemFactory(final String ctype, final CoordinateReferenceFrame refSystem) {
        final AbstractCrs crs;
        switch (ctype) {
            case "RA":
            case "DEC":
                crs = createEquatorial(refSystem);
                break;
            case "GLON":
            case "GLAT":
                crs = new Galactic();
                break;
            case "ELON":
            case "ELAT":
                crs = createEcliptic(refSystem);
                break;
            default:
                throw new JWcsError("The coordinate reference system (" + ctype + "," + refSystem + ") is not supported");
        }
        return crs;
    }

    /**
     * Creates the equatorial coordinate reference system.
     *
     * @param refSystem the coordinate reference frame
     * @return the equatorial coordinate reference system
     */
    private AbstractCrs createEquatorial(final CoordinateReferenceFrame refSystem) {
        final AbstractCrs crs = new Equatorial();
        if (refSystem != null) {
            crs.setCoordinateReferenceFrame(refSystem);
        }
        return crs;
    }

    /**
     * Create the ecliptic coordinate reference system.
     *
     * @param refSystem the coordinate reference frame
     * @return the ecliptic coordinate reference system
     */
    private AbstractCrs createEcliptic(final CoordinateReferenceFrame refSystem) {
        final AbstractCrs crs = new Ecliptic();
        if (refSystem != null) {
            crs.setCoordinateReferenceFrame(refSystem);
        }
        return crs;
    }

    /**
     * Make the initialization of the WCS. 
     * 
     * <p>By default, you must call
     * {@link #init()}
     *
     * @throws io.github.malapert.jwcs.proj.exception.JWcsException when an
     * error occurs
     */
    public abstract void doInit() throws JWcsException;

    /**
     * Returns the parameters of the projection.
     *
     * @return the projection parameters.
     */
    public final ProjectionParameter[] getProjectionParameters() {
        return getProj().getProjectionParameters();
    }

    /**
     * Returns the projection's name.
     *
     * @return the projection's name
     */
    public final String getName() {
        return getProj().getName();
    }

    /**
     * Returns the projection family.
     * 
     * <p>The supported projection families are the following:
     * <ul>
     * <li>{@link io.github.malapert.jwcs.proj.AbstractCylindricalProjection}</li>
     * <li>{@link io.github.malapert.jwcs.proj.AbstractConicProjection}</li>
     * <li>{@link io.github.malapert.jwcs.proj.AbstractPolyConicProjection}</li>
     * <li>{@link io.github.malapert.jwcs.proj.AbstractZenithalProjection}</li>
     * </ul>
     *
     * @return the projection's name
     */
    public final String getNameFamily() {
        return getProj().getNameFamily();
    }

    /**
     * Returns the projection's description.
     *
     * @return the projection's description
     */
    public final String getDescription() {
        return getProj().getDescription();
    }

    /**
     * Computes CD matrix from CDELT[] and CROTA.
     *
     * <p>The computation is realized as follows:
     * <pre>
     *   cos0 = cos(crota)
     *   sin0 = sin(crota)
     *   cd11 = cdelt1 * cos0
     *   cd12 = abs(cdelt2) * signum(cdelt1) * sin0
     *   cd21 = -abs(cdelt1) * signum(cdelt2) * sin0;
     *   cd22 = cdelt2 * cos0;
     * </pre>
     *
     * @param cdelt increment of position
     * @param crota rotation
     * @return the cd matrix as array
     */
    protected static final double[][] computeCdFromCdelt(final double[] cdelt, final double crota) {
        final double cos0 = FastMath.cos(FastMath.toRadians(crota));
        final double sin0 = FastMath.sin(FastMath.toRadians(crota));
        final double cd11 = cdelt[0] * cos0;
        final double cd12 = FastMath.abs(cdelt[1]) * FastMath.signum(cdelt[0]) * sin0;
        final double cd21 = -FastMath.abs(cdelt[0]) * FastMath.signum(cdelt[1]) * sin0;
        final double cd22 = cdelt[1] * cos0;
        final double[][] array = {
            {cd11, cd12},
            {cd21, cd22}
        };
        return array;
    }

    /**
     * Computes CD matrix from PC[][] and CDELT[].
     *
     * @param pc pc matrix
     * @param cdelt cdelt array
     * @return the CD matrix
     */
    protected static double[][] pc2cd(final double[][] pc, final double[] cdelt) {
        final double[][] cd_conv = {
            {cdelt[0] * pc[0][0], cdelt[1] * pc[1][0]},
            {cdelt[0] * pc[0][1], cdelt[1] * pc[1][1]}
        };
        return cd_conv;
    }

    /**
     * Creates the CD matrix.
     * 
     * <p>The CD matrix is created by reading CD matrix or by computing the CD
     * matrix from the CDELT and CROTA.
     *
     * @return the CD matrix
     */
    protected final RealMatrix createCdMatrix() {
        final double[][] arraycd;
        if (hasCd()) {
            arraycd = new double[2][2];
            arraycd[0][0] = cd(1, 1);
            arraycd[0][1] = cd(1, 2);
            arraycd[1][0] = cd(2, 1);
            arraycd[1][1] = cd(2, 2);
        } else {
            final double[] cdelt = new double[]{getValueAsDouble(CDELT1), getValueAsDouble(CDELT2)};
            arraycd = computeCdFromCdelt(cdelt, getValueAsDouble(CROTA2));
        }
        return createRealMatrix(arraycd);
    }

    @Override
    public boolean hasCd() {
        return hasKeyword(CD11)
                || (hasKeyword(CDELT1) && hasKeyword(CROTA2))
                || (hasKeyword(CDELT1) && hasKeyword(PC11));
    }

    @Override
    public double cd(final int i, final int j) {
        final double result;
        if (hasKeyword(CD11)) {
            result = this.getValueAsDouble("CD" + i + "_" + j);
        } else if (hasKeyword(CROTA2)) {
            final double[] cdelt = new double[]{getValueAsDouble(CDELT1), getValueAsDouble(CDELT2)};
            final double[][] cdTmp = AbstractJWcs.computeCdFromCdelt(cdelt, getValueAsDouble(CROTA2));
            result = cdTmp[i - 1][j - 1];
        } else if (hasKeyword(PC11)) {
            final double[][] pc = new double[][]{{getValueAsDouble(PC11), getValueAsDouble(PC12)},
            {getValueAsDouble(PC21), getValueAsDouble(PC22)}};
            final double[] cdelt = new double[]{getValueAsDouble(CDELT1), getValueAsDouble(CDELT2)};
            final double[][] cdTmp = AbstractJWcs.pc2cd(pc, cdelt);
            result = cdTmp[i - 1][j - 1];
        } else {
            throw new JWcsError("cd" + i + j + " not found");
        }
        return result;
    }

    @Override
    public double lonpole() {
        final double result;
        if (hasKeyword(LONPOLE)) {
            result = getValueAsDouble(LONPOLE);
        } else if (hasKeyword(PV13)) {
            result = getValueAsDouble(PV13);
        } else {
            result = Double.NaN;
        }
        return result;
    }

    @Override
    public double latpole() {
        final double result;
        if (hasKeyword(LATPOLE)) {
            result = getValueAsDouble(LATPOLE);
        } else if (hasKeyword(PV14)) {
            result = getValueAsDouble(PV14);
        } else {
            result = Double.NaN;
        }
        return result;
    }

    @Override
    public int wcsaxes() {
        if (hasKeyword(NAXIS)) {
            return getValueAsInt(NAXIS);
        } else {
            throw new JWcsError(NAXIS + " not found");
        }
    }

    @Override
    public int naxis(final int j) {
        if (hasKeyword("NAXIS" + j)) {
            return getValueAsInt("NAXIS" + j);
        } else {
            throw new JWcsError("NAXIS" + j + " not found");
        }
    }

    @Override
    public double crval(final int n) {
        if (hasKeyword("CRVAL" + n)) {
            return getValueAsDouble("CRVAL" + n);
        } else {
            throw new JWcsError("CRVAL" + n + " not found");
        }
    }

    @Override
    public double crpix(final int n) {
        if (hasKeyword("CRPIX" + n)) {
            return getValueAsDouble("CRPIX" + n);
        } else {
            throw new JWcsError("CRPIX" + n + " not found");
        }
    }

    @Override
    public String ctype(final int n) {
        if (hasKeyword("CTYPE" + n)) {
            return getValueAsString("CTYPE" + n);
        } else {
            throw new JWcsError("CTYPE" + n + " not found");
        }
    }

    @Override
    public double pv(final int i, final int m) {
        return getValueAsDouble("PV" + i + "_" + m);
    }

    @Override
    public String cunit(final int i) {
        return getValueAsString("CUNIT" + i);
    }

    @Override
    public double equinox() {
        return getValueAsDouble(EQUINOX);
    }

    /**
     * Returns the value of keyword when keyword exists otherwise defaultValue.
     *
     * @param keyword the keyword
     * @param defaultValue the default value
     * @return the value
     */
    public double getValueAsDouble(final String keyword, final double defaultValue) {
        final double result;
        if (hasKeyword(keyword)) {
            result = getValueAsDouble(keyword);
        } else {
            LOG.log(Level.WARNING, "{0} not found -- use default value {1}", new Object[]{keyword, defaultValue});
            result = defaultValue;
        }
        return result;
    }

    /**
     * Scale factor. 
     * 
     * <p>This is used to convert the CRVAL into degree.
     *
     * @param cunit The cunit axis
     * @return the scale factor to apply at CRVAL
     */
    private double convertToDegree(final String cunit) {
        final double cx;
        if (hasKeyword(cunit)) {
            final String unit_lc = cunit.toLowerCase(Locale.ENGLISH);
            switch (unit_lc) {
                case "acmin":
                    cx = 1 / 60.0;
                    break;
                case "arcsec":
                    cx = 1 / 3600.0;
                    break;
                case "mas":
                    cx = 1 / 3600000.0;
                    break;
                case "rad":
                    cx = 180 / FastMath.PI;
                    break;
                default:
                    cx = 1;
                    break;
            }
        } else {
            // assumes it is in degree;
            cx = 1;
        }
        return cx;
    }

    /**
     * Creates the projection by reading CTYPE1.
     * 
     * @return the projection
     * @throws BadProjectionParameterException when the projection parameter is
     * wrong
     */
    protected final AbstractProjection createProjection() throws BadProjectionParameterException {
        final String ctype1 = ctype(1);
        final String codeProjection = ctype1.substring(ctype1.lastIndexOf('-') + 1, ctype1.length());
        final double cx = convertToDegree(cunit(1));
        final double cy = convertToDegree(cunit(2));
        final AbstractProjection projection = createProjectionFactory(codeProjection, cx, cy);
        setNativeLongitudeOfFiducialPoint(projection);
        setNativeLatitudeOfFiducialPoint(projection);
        setNativeLongitudeOfCelestialPole(projection);
        setNativeLatitudeOfCelestialPole(projection);
        //TO DO  apply shift PV10
        return projection;
    }

    /**
     * Creates a projection based on its projection code.
     *
     * @param projectionCode projection code
     * @param cx scale factor along X
     * @param cy scale factor along Y
     * @return the projection
     * @throws BadProjectionParameterException When a bad parameter is provided
     * to the projection
     * @throws JWcsError projection code is not supported
     */
    private AbstractProjection createProjectionFactory(final String projectionCode, final double cx, final double cy) throws BadProjectionParameterException {
        final AbstractProjection projection;
        switch (projectionCode) {
            case "ZPN":
                projection = createZPNProjection(cx, cy);
                break;
            case "TPV":
                projection = createTPVProjection(cx, cy);
                break;
            case "BON":
                LOG.log(Level.INFO, "Creates a AIT projection with (crval1,crval2)=({0},{1}) theta1={2}", new Object[]{crval(1) * cx, crval(2) * cx, getValueAsDouble(PV21, 0)});
                projection = new BON(crval(1) * cx, crval(2) * cy, getValueAsDouble(PV21, 0));
                break;
            case "CEA":
                LOG.log(Level.INFO, "Creates a CEA projection with (crval1,crval2)=({0},{1}) lambda={2}", new Object[]{crval(1) * cx, crval(2) * cx, getValueAsDouble(PV21, 0)});
                projection = new CEA(crval(1) * cx, crval(2) * cy, getValueAsDouble(PV21, 1));
                break;                
            case "COD":
                LOG.log(Level.INFO, "Creates a COD projection with (crval1,crval2)=({0},{1}) (theta_a,eta)=({2},{3})", new Object[]{crval(1) * cx, crval(2) * cx, getValueAsDouble(PV21, 0), getValueAsDouble(PV22, 0)});
                projection = new COD(crval(1) * cx, crval(2) * cy, getValueAsDouble(PV21, 0), getValueAsDouble(PV22, 0));
                break;
            case "COE":
                LOG.log(Level.INFO, "Creates a COE projection with (crval1,crval2)=({0},{1}) (theta_a,eta)=({2},{3})", new Object[]{crval(1) * cx, crval(2) * cx, getValueAsDouble(PV21, 0), getValueAsDouble(PV22, 0)});
                projection = new COE(crval(1) * cx, crval(2) * cy, getValueAsDouble(PV21, 0), getValueAsDouble(PV22, 0));
                break;
            case "COO":
                LOG.log(Level.INFO, "Creates a COO projection with (crval1,crval2)=({0},{1}) (theta_a,eta)=({2},{3})", new Object[]{crval(1) * cx, crval(2) * cx, getValueAsDouble(PV21, 0), getValueAsDouble(PV22, 0)});
                projection = new COO(crval(1) * cx, crval(2) * cy, getValueAsDouble(PV21, 0), getValueAsDouble(PV22, 0));
                break;
            case "COP":
                LOG.log(Level.INFO, "Creates a COP projection with (crval1,crval2)=({0},{1}) (theta_a,eta)=({2},{3})", new Object[]{crval(1) * cx, crval(2) * cx, getValueAsDouble(PV21, 0), getValueAsDouble(PV22, 0)});
                projection = new COP(crval(1) * cx, crval(2) * cy, getValueAsDouble(PV21, 0), getValueAsDouble(PV22, 0));
                break;                               
            case "SZP":
                projection = createSZPProjection(cx, cy);
                break;
            case "NCP":
                LOG.log(Level.INFO, "Creates a NCP projection with (crval1,crval2)=({0},{1}) (ksi,eta)=({2},{3})", new Object[]{crval(1) * cx, crval(2) * cx, 0, 1/FastMath.tan(getValueAsDouble(CRVAL2))});
                projection = new SIN(crval(1) * cx, crval(2) * cy, 0, 1/FastMath.tan(getValueAsDouble(CRVAL2)));
                break;
            default:
                projection = createProjection(projectionCode, cx, cy);
                break;
        }
        return projection;
    }

    /**
     * Creates dynamically a projection based on the projection code.
     * 
     * <p>No projection parameter is used to instantiate this class.
     *
     * @param projectionCode projection to instantiate dynamically
     * @param cx scale factor along X
     * @param cy scale factor along Y
     * @return the projection corresponding to the projection code
     * @throws JWcsError Cannot find the projection or error when instantiate
     * dynamically the projection
     */
    private AbstractProjection createStandardProjection(final String projectionCode, final double cx, final double cy) {
        final AbstractProjection projection;
        try {
            final Package packageName = this.getClass().getPackage();
            final String name = packageName.getName()+".proj.";
            final Class<?> clazz = Class.forName(name + projectionCode);
            final Constructor<?> constructor = clazz.getConstructor(Double.TYPE, Double.TYPE);
            final Object instance = constructor.newInstance(crval(1) * cx, crval(2) * cy);            
            projection = (AbstractProjection) instance;
            LOG.log(Level.INFO, "Creates a {0} projection with (crval1,crval2)=({1},{2})", new Object[]{projectionCode, crval(1) * cx, crval(2) * cx});
        } catch (ClassNotFoundException ex) {
            throw new JWcsError("The projection " + projectionCode + " is not supported.");
        } catch (NoSuchMethodException | SecurityException | InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException ex) {
            throw new JWcsError(ex);
        }
        return projection;
    }

    /**
     * Creates dynamically a projection based on the projection code.
     * 
     * <p>projection parameters are used to instantiate this class.
     *
     * @param projectionCode projection to instantiate dynamically
     * @param cx scale factor along X
     * @param cy scale factor along Y
     * @return the projection corresponding to the projection code
     * @throws JWcsError Cannot find the projection or error when instantiate
     * dynamically the projection
     */
    private AbstractProjection createStandardProjectionWithParameters(final String projectionCode, final double cx, final double cy) {
        final AbstractProjection projection;
        final Package packageName = this.getClass().getPackage();
        final String name = packageName.getName() + ".proj.";
        final double pv21 = getValueAsDouble(PV21);        
        try {
            final Class<?> clazz = Class.forName(name + projectionCode);
            final Constructor<?> constructor;
            final Object instance;
            if (hasKeyword(PV22)){
                final double pv22 = getValueAsDouble(PV22);
                constructor = clazz.getConstructor(Double.TYPE, Double.TYPE, Double.TYPE, Double.TYPE);                
                instance = constructor.newInstance(crval(1) * cx, crval(2) * cy, pv21, pv22);
                LOG.log(Level.INFO, "Creates a {0} projection with (crval1,crval2)=({1},{2}) (pv21,pv22)=({3},{4})", new Object[]{projectionCode, crval(1) * cx, crval(2) * cx, pv21, pv22});                
            } else {
                constructor = clazz.getConstructor(Double.TYPE, Double.TYPE, Double.TYPE);
                instance = constructor.newInstance(crval(1) * cx, crval(2) * cy, pv21);
                LOG.log(Level.INFO, "Creates a {0} projection with (crval1,crval2)=({1},{2}) pv21={3}", new Object[]{projectionCode, crval(1) * cx, crval(2) * cx, pv21});
            }
            projection = (AbstractProjection) instance;
        } catch (ClassNotFoundException ex) {
            throw new JWcsError("The projection " + projectionCode + " is not supported.");
        } catch (NoSuchMethodException | SecurityException | InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException ex) {
            throw new JWcsError(ex);
        }
        return projection;
    }

    /**
     * Creates a projection.
     *
     * <p>Instantiates either projection with projection parameters when they are
     * available otherwise instantiates only the projection with no projection 
     * parameter.
     * 
     * @param projectionCode projection to instantiate dynamically
     * @param cx scale factor along X
     * @param cy scale factor along Y
     * @return the projection corresponding to the projection code
     * @throws JWcsError Cannot find the projection or error when instantiate
     * dynamically the projection
     */
    private AbstractProjection createProjection(final String projectionCode, final double cx, final double cy) {
        final AbstractProjection projection;    
         if (hasKeyword(PV21) && hasKeyword(PV22)) {
             projection = createStandardProjectionWithParameters(projectionCode, cx, cy);
         } else if (hasKeyword(PV21)) {
             projection = createStandardProjectionWithParameters(projectionCode, cx, cy);
         } else {
             projection = createStandardProjection(projectionCode, cx, cy);
         }
         return projection;
    }
    
    /**
     * Creates a SZP projection.
     *
     * @param cx scale factor along X
     * @param cy scale factor along Y
     * @return the SZP projection
     * @throws BadProjectionParameterException when a bad parameter is provided
     * to the projection
     */
    private AbstractProjection createSZPProjection(final double cx, final double cy) throws BadProjectionParameterException {
        final AbstractProjection projection;
        if (hasKeyword(PV21) && hasKeyword(PV22) && hasKeyword(PV23)) {
            LOG.log(Level.INFO, "Creates a SZP projection with (crval1,crval2)=({0},{1}) (mu,phic,thetac)=({2},{3},{4})", new Object[]{crval(1) * cx, crval(2) * cx, getValueAsDouble(PV21), getValueAsDouble(PV22), getValueAsDouble(PV23)});
            projection = new SZP(crval(1) * cx, crval(2) * cy, getValueAsDouble(PV21), getValueAsDouble(PV22), getValueAsDouble(PV23));
        } else {
            projection = createStandardProjection("SZP", crval(1) * cx, crval(2) * cy);
        }
        return projection;
    }

    /**
     * Creates ZPN projection.
     *
     * @param cx the scale factor along X
     * @param cy the scale factor along Y
     * @return the ZPN projection
     * @throws BadProjectionParameterException when a bad parameter is provided
     * to the projection
     */
    private AbstractProjection createZPNProjection(final double cx, final double cy) throws BadProjectionParameterException {
        final Iterator iter = iterator();
        final Map<String, Double> pvMap = new HashMap<>();
        while (iter.hasNext()) {
            final Object keyObj = iter.next();
            if (keyObj instanceof HeaderCard) {
                final HeaderCard card = (HeaderCard) keyObj;
                final String key = card.getKey();
                if (key.startsWith("PV2")) {
                    pvMap.put(key, getValueAsDouble(key));
                }
            } else {
                final String key = (String) keyObj;
                if (key.startsWith("PV2")) {
                    pvMap.put(key, getValueAsDouble(key));
                }
            }
        }
        final double[] pvsPrimitif = new double[pvMap.size()];
        for (int i = 0; i < pvMap.size(); i++) {
            pvsPrimitif[i] = pvMap.get("PV2_" + i);
        }
        LOG.log(Level.INFO, "Creates a ZPN projection with (crval1,crval2)=({0},{1} PV={2})", new Object[]{crval(1) * cx, crval(2) * cx, Arrays.toString(pvsPrimitif)});
        return new ZPN(crval(1) * cx, crval(2) * cy, pvsPrimitif);
    }

    /**
     * Creates TPV projection.
     *
     * @param cx the scale factor along X
     * @param cy the scale factor along Y
     * @return the TPV projection
     * @throws BadProjectionParameterException when a bad parameter is provided
     * to the projection
     */
    private AbstractProjection createTPVProjection(final double cx, final double cy) throws BadProjectionParameterException {
        final Iterator iter = iterator();
        final Map<String, Double> pvMap_a = new HashMap<>();
        final Map<String, Double> pvMap_b = new HashMap<>();
        while (iter.hasNext()) {
            final Object keyObj = iter.next();
            if (keyObj instanceof HeaderCard) {
                final HeaderCard card = (HeaderCard) keyObj;
                final String key = card.getKey();
                if (key.startsWith("PV1_"))
                    pvMap_a.put(key, getValueAsDouble(key));
                if (key.startsWith("PV2_"))
                    pvMap_b.put(key, getValueAsDouble(key));
            } else {
                final String key = (String) keyObj;
                if (key.startsWith("PV1_"))
                    pvMap_a.put(key, getValueAsDouble(key));
                if (key.startsWith("PV2_"))
                    pvMap_b.put(key, getValueAsDouble(key));
            }
        }
        
        int pv_size = 40;
        final double[] pv_a = new double[pv_size];
        for (int i = 0; i < pv_a.length; i++) {
            pv_a[i] = (i == 1) ? 1 : 0;
            String key = "PV1_" + i;
            if (pvMap_a.containsKey(key))
                pv_a[i] = pvMap_a.get(key);
        }
        final double[] pv_b = new double[pv_size];
        for (int i = 0; i < pv_b.length; i++) {
            pv_b[i] = (i == 1) ? 1 : 0;
            String key = "PV2_" + i;
            if (pvMap_b.containsKey(key))
                pv_b[i] = pvMap_b.get(key);
        }
        
        LOG.log(Level.INFO, "Creates a TPV projection with (crval1,crval2)=({0},{1} PV_a={2} PV_b={3})", new Object[]{crval(1) * cx, crval(2) * cx, Arrays.toString(pv_a), Arrays.toString(pv_b)});
        return new TPV(crval(1) * cx, crval(2) * cy, pv_a, pv_b);
    }
    
    /**
     * Sets the Native longitude of the fiducial point to the projection.
     *
     * @param projection the projection
     */
    private void setNativeLongitudeOfFiducialPoint(final AbstractProjection projection) {
        if (hasKeyword(PV11)) {
            LOG.log(Level.INFO, "Sets phi0 to {0}", getValueAsDouble(PV11));
            projection.setPhi0(getValueAsDouble(PV11));
        }
    }

    /**
     * Sets the Native latitude of the celestial pole to the projection.
     *
     * @param projection the projection
     */
    private void setNativeLatitudeOfFiducialPoint(final AbstractProjection projection) {
        if (hasKeyword(PV12)) {
            LOG.log(Level.INFO, "Sets theta0 to {0}", getValueAsDouble(PV12));
            projection.setTheta0(getValueAsDouble(PV12));
        }
    }

    /**
     * Sets the Native longitude of the celestial pole to the projection.
     *
     * @param projection the projection
     */
    private void setNativeLongitudeOfCelestialPole(final AbstractProjection projection) {
        if (!Double.isNaN(lonpole())) {
            LOG.log(Level.INFO, "Sets phip to {0}", lonpole());
            projection.setPhip(FastMath.toRadians(lonpole()));
        }
    }

    /**
     * Sets the Native latitude of the celestial pole to the projection.
     *
     * @param projection the projection
     */
    private void setNativeLatitudeOfCelestialPole(final AbstractProjection projection) {
        if (!Double.isNaN(latpole())) {
            LOG.log(Level.INFO, "Sets thetap to {0}", latpole());
            projection.setThetap(FastMath.toRadians(latpole()));
        }
    }

    /**
     * Transforms the position of a pixel given by (x,y) in a position in the
     * sky.
     *
     * @param x X coordinate of the pixel
     * @param y Y coordinate of the pixel
     * @return the pixel position in the sky
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException when
     * there is a projection error
     */
    @Override
    public double[] pix2wcs(final double x, final double y) throws ProjectionException {
        final double[][] arraypj = {
            {x - crpix(1), y - crpix(2)}
        };
        final RealMatrix pj = createRealMatrix(arraypj);
        final RealMatrix pjt = pj.transpose();
        final RealMatrix xi = this.getCd().multiply(pjt);
        return this.getProj().projectionPlane2wcs(xi.getEntry(0, 0), xi.getEntry(1, 0));
    }

    /**
     * Transforms an array of pixel position in an array of position in the sky.
     *
     * @param pixels an array of pixel
     * @return an array of sky position
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException when
     * there is a projection error
     * @throws JWcsError the length of pixels must be a multiple of 2
     */
    @Override
    public double[] pix2wcs(final double[] pixels) throws ProjectionException {
        final int pixelsLength = pixels.length;
        if (pixelsLength % 2 != 0) {
            throw new JWcsError("the length of pixels must be a multiple of 2");
        }
        final double[] skyPositions = new double[pixelsLength];
        for (int i = 0; i < pixelsLength; i = i + 2) {
            final double x = pixels[i];
            final double y = pixels[i + 1];
            final double[] result = this.pix2wcs(x, y);
            skyPositions[i] = result[0];
            skyPositions[i + 1] = result[1];
        }
        return skyPositions;
    }

    @Override
    public double[] getCenter() throws ProjectionException {
        return pix2wcs(0.5 * naxis(1), 0.5 * naxis(2));
    }

    @Override
    public double[] getFov() throws ProjectionException {
        return pix2wcs(new double[]{0.5, 0.5, naxis(1) + 0.5, 0.5, naxis(1) + 0.5, naxis(2) + 0.5, 0.5, naxis(2) + 0.5});
    }

    /**
     * Checks validity of longitude and latitude.
     *
     * @param longitude longitude [0, 360]
     * @param latitude latitude [-90, 90]
     * @throws JWcsError the range is not valid
     */
    private void checkLongitudeLatitude(final double longitude, final double latitude) {
        if (longitude > MAX_LONGITUDE || longitude < MIN_LONGITUDE) {
            throw new JWcsError("Longitude must be [0, 360], found " + longitude);
        }
        if (latitude > MAX_LATITUDE || latitude < MIN_LATITUDE) {
            throw new JWcsError("Latitude must be [-90, 90], found " + latitude);
        }
    }

    /**
     * Returns true if the given lat/lon point is visible in this projection.
     *
     * @param lon longitude in degrees.
     * @param lat latitude in degrees.
     * @return True when the point is visible otherwise False.
     */
    @Override
    public boolean inside(final double lon, final double lat) {
        return this.getProj().inside(FastMath.toRadians(lon), FastMath.toRadians(lat));
    }

    /**
     * Checks if the line is visible.
     *
     * @param pos1 first point of the line
     * @param pos2 last point of the line
     * @return True when the line is visible otherwise False.
     */
    public boolean isLineToDraw(final double[] pos1, final double[] pos2) {
        final boolean result;
        final boolean isFinite = Double.isFinite(pos1[0]) && Double.isFinite(pos1[1]) && Double.isFinite(pos2[0]) && Double.isFinite(pos2[1]);
        if (isFinite) {
            result = this.getProj().isLineToDraw(pos1, pos2);
        } else {
            result = false;
        }
        return result;
    }

    /**
     * Transforms the sky position given by (longitude, latitude) in a pixel
     * position.
     *
     * @param longitude longitude of the sky position
     * @param latitude latitude of the sky position
     * @return the sky position in the pixel grid.
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException when
     * there is a projection error
     */
    @Override
    public double[] wcs2pix(final double longitude, final double latitude) throws ProjectionException {
        checkLongitudeLatitude(longitude, latitude);
        final double[] coordVal = this.getProj().wcs2projectionPlane(FastMath.toRadians(longitude), FastMath.toRadians(latitude));
        final double[][] coord = {
            {coordVal[0], coordVal[1]}
        };
        final RealMatrix coordM = createRealMatrix(coord);
        final RealMatrix matrix = coordM.multiply(getCdInverse());
        return new double[]{matrix.getEntry(0, 0) + crpix(1), matrix.getEntry(0, 1) + crpix(2)};
    }

    /**
     * Transforms an array of sky position in an array of pixel position.
     *
     * @param skyPositions array of sky positions
     * @return the sky position in the pixel grid.
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException when 
     * there is a projection error
     * @throws JWcsError When the length of <code>skyPositions</code> is not a
     * multiple of 2
     */
    @Override
    public double[] wcs2pix(final double[] skyPositions) throws ProjectionException {
        final int skyPositionLength = skyPositions.length;
        if (skyPositionLength % 2 != 0) {
            throw new JWcsError("the length of skyPositions must be a multiple of 2");
        }
        final double[] pixelPositions = new double[skyPositionLength];
        for (int i = 0; i < skyPositionLength; i = i + 2) {
            final double ra = skyPositions[i];
            final double dec = skyPositions[i + 1];
            final double[] result = this.wcs2pix(ra, dec);
            pixelPositions[i] = result[0];
            pixelPositions[i + 1] = result[1];
        }
        return pixelPositions;
    }

    /**
     * Returns the projection.
     *
     * @return the projection
     */
    protected final AbstractProjection getProj() {
        return proj;
    }

    /**
     * Sets the projection.
     *
     * @param proj the projection to set
     */
    protected final void setProj(final AbstractProjection proj) {
        this.proj = proj;
    }

    /**
     * Returns the CD matrix.
     *
     * @return the cd
     */
    protected final RealMatrix getCd() {
        return cd;
    }

    /**
     * The CD matrix to set.
     *
     * @param cd the cd to set
     */
    protected final void setCd(final RealMatrix cd) {
        this.cd = cd;
    }

    /**
     * Returns the CD matrix inverse.
     *
     * @return the CD matrix inverse
     */
    protected final RealMatrix getCdInverse() {
        return cdInverse;
    }

    /**
     * Sets the CD matrix inverse.
     *
     * @param cdInverse the CD matrix inverse to set
     */
    protected final void setCdInverse(final RealMatrix cdInverse) {
        this.cdInverse = cdInverse;
    }
}
