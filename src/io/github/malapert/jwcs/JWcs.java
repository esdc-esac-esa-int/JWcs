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
package io.github.malapert.jwcs;

import io.github.malapert.jwcs.coordsystem.Ecliptic;
import io.github.malapert.jwcs.coordsystem.Equatorial;
import io.github.malapert.jwcs.coordsystem.FK4;
import io.github.malapert.jwcs.coordsystem.FK4_NO_E;
import io.github.malapert.jwcs.coordsystem.FK5;
import io.github.malapert.jwcs.coordsystem.Galactic;
import io.github.malapert.jwcs.coordsystem.ICRS;
import io.github.malapert.jwcs.coordsystem.Crs;
import io.github.malapert.jwcs.proj.AIT;
import io.github.malapert.jwcs.proj.ARC;
import io.github.malapert.jwcs.proj.AZP;
import io.github.malapert.jwcs.proj.BON;
import io.github.malapert.jwcs.proj.CAR;
import io.github.malapert.jwcs.proj.CEA;
import io.github.malapert.jwcs.proj.COD;
import io.github.malapert.jwcs.proj.COE;
import io.github.malapert.jwcs.proj.COO;
import io.github.malapert.jwcs.proj.COP;
import io.github.malapert.jwcs.proj.CYP;
import io.github.malapert.jwcs.proj.MER;
import io.github.malapert.jwcs.proj.MOL;
import io.github.malapert.jwcs.proj.PAR;
import io.github.malapert.jwcs.proj.PCO;
import io.github.malapert.jwcs.proj.Projection;
import io.github.malapert.jwcs.proj.Projection.ProjectionParameter;
import io.github.malapert.jwcs.proj.SFL;
import io.github.malapert.jwcs.proj.SIN;
import io.github.malapert.jwcs.proj.STG;
import io.github.malapert.jwcs.proj.SZP;
import io.github.malapert.jwcs.proj.TAN;
import io.github.malapert.jwcs.proj.ZEA;
import io.github.malapert.jwcs.proj.ZPN;
import io.github.malapert.jwcs.proj.exception.BadProjectionParameterException;
import io.github.malapert.jwcs.proj.exception.JWcsException;
import io.github.malapert.jwcs.proj.exception.JWcsError;
import io.github.malapert.jwcs.proj.exception.ProjectionException;
import io.github.malapert.jwcs.utility.TimeUtils;
import java.text.ParseException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import nom.tam.fits.Header;
import nom.tam.fits.HeaderCard;
import io.github.malapert.jwcs.coordsystem.CoordinateReferenceFrame;
import static io.github.malapert.jwcs.utility.NumericalUtils.createRealMatrix;
import static io.github.malapert.jwcs.utility.NumericalUtils.inverse;
import org.apache.commons.math3.linear.RealMatrix;

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
public abstract class JWcs implements JWcsKeyProvider {

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
     * Number of axes. 2 for an image
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
     * Reference along X axis in pixel frame. This keyword is required for
     * projection computation.
     */
    public static final String CRPIX1 = "CRPIX1";
    /**
     * Reference along Y axis in pixel frame. This keyword is required for
     * projection computation.
     */
    public static final String CRPIX2 = "CRPIX2";
    /**
     * Reference along longitude in degrees in celestial frame. This keyword is
     * required for projection computation.
     */
    public static final String CRVAL1 = "CRVAL1";
    /**
     * Reference along latitude in degrees in celestial frame. This keyword is
     * required for projection computation.
     */
    public static final String CRVAL2 = "CRVAL2";
    /**
     * Projection type along X axis. This keyword is required for projection
     * computation.
     */
    public static final String CTYPE1 = "CTYPE1";
    /**
     * Projection type along Y axis. This keyword is required for projection
     * computation.
     */
    protected static final String CTYPE2 = "CTYPE2";
    /**
     * Scale (degrees / pixel) and rotation matrix. For projection computation,
     * information about scale and rotation are needed. Either the CD matrix is
     * provided or the following element (CDELT1, CDELT2, CROTA2) or (PC matrix,
     * CDELT1, CDELT2).
     */
    public static final String CD11 = "CD1_1";
    /**
     * Scale (degrees / pixel) and rotation matrix. For projection computation,
     * information about scale and rotation are needed. Either the CD matrix is
     * provided or the following element (CDELT1, CDELT2, CROTA2) or (PC matrix,
     * CDELT1, CDELT2).
     */
    public static final String CD12 = "CD1_2";
    /**
     * Scale (degrees / pixel) and rotation matrix. For projection computation,
     * information about scale and rotation are needed. Either the CD matrix is
     * provided or the following element (CDELT1, CDELT2, CROTA2) or (PC matrix,
     * CDELT1, CDELT2).
     */
    public static final String CD21 = "CD2_1";
    /**
     * Scale (degrees / pixel) and rotation matrix. For projection computation,
     * information about scale and rotation are needed. Either the CD matrix is
     * provided or the following element (CDELT1, CDELT2, CROTA2) or (PC matrix,
     * CDELT1, CDELT2).
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
     * Scale (degrees / pixel) along X axis when CD matrix is not defined. For
     * projection computation, information about scale and rotation are needed.
     * Either the CD matrix is provided or the following element (CDELT1,
     * CDELT2, CROTA2) or (PC matrix, CDELT1, CDELT2).
     */
    public static final String CDELT1 = "CDELT1";
    /**
     * Scale (degrees / pixel) along X axis when CD matrix is not defined. For
     * projection computation, information about scale and rotation are needed.
     * Either the CD matrix is provided or the following element (CDELT1,
     * CDELT2, CROTA2) or (PC matrix, CDELT1, CDELT2).
     */
    public static final String CDELT2 = "CDELT2";
    /**
     * For projection computation, information about scale and rotation are
     * needed. Either the CD matrix is provided or the following element
     * (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1, CDELT2).
     */
    public static final String CROTA2 = "CROTA2";
    /**
     * Equinox value.
     */
    public static final String EQUINOX = "EQUINOX";
    /**
     * Deformation matrix. For projection computation, information about scale
     * and rotation are needed. Either the CD matrix is provided or the
     * following element (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1,
     * CDELT2).
     */
    public static final String PC11 = "PC1_1";
    /**
     * Deformation matrix. For projection computation, information about scale
     * and rotation are needed. Either the CD matrix is provided or the
     * following element (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1,
     * CDELT2).
     */
    public static final String PC12 = "PC1_2";
    /**
     * Deformation matrix. For projection computation, information about scale
     * and rotation are needed. Either the CD matrix is provided or the
     * following element (CDELT1, CDELT2, CROTA2) or (PC matrix, CDELT1,
     * CDELT2).
     */
    public static final String PC21 = "PC2_1";
    /**
     * Deformation matrix. For projection computation, information about scale
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
     * Projection object.
     */
    private Projection proj;

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
    protected static final Logger LOG = Logger.getLogger(JWcs.class.getName());

    /**
     * Initialize the WCS Object.
     *
     * <p>
     * The WCS object is initialized by doing the following steps:
     * <ul>
     * <li>creates the projection</li>
     * <li>creates the CD matrix</li>
     * <li>creates the CD matrix inverse</li>
     * <li>checks the WCS</li>
     * </ul>
     *
     * @throws io.github.malapert.jwcs.proj.exception.JWcsException When WCS is
     * not valid
     */
    protected final void init() throws JWcsException {
        checkWcs();
        setProj(createProjection());
        setCd(createCdMatrix());
        setCdInverse(inverse(getCd()));
    }

    /**
     * Checks WCS keywords.
     * <p>
     * Raises an exception when a problem is detected.
     * </p>
     *
     * @throws JWcsException When WCS is not valid
     */
    protected void checkWcs() throws JWcsException {
        // do nothing
    }

    /**
     * Checks if the WCS header is valid.
     *
     * @param hdr Header FITS
     * @return True when the Header is valid otherwise False.
     */
    public static boolean isValidWcs(final Header hdr) {
        final JWcs wcs = new JWcsFits(hdr);
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
     * <p>
     * Returns the value of the MJD-OBS keyword when it is present otherwise
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
                mjd = String.valueOf(TimeUtils.convertISOToModifiedJulianDate(getValueAsString("DATE-OBS")));
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
     * <p>
     * To find the reference system, the algorithm proceed as it:
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
     * Create a FK4 coordinate reference frame with a Besselian epoch of equinox
     * and a Modified Julian date as epoch of observation when equinox is smaller
     * then 1984. Otherwise, a FK5 coordinate refrence frame is created based on
     * a Julian date as epoch of observation.
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
     * Creates a factory to create the right coordinate reference frame based
     * on the RADESYS keyword value.
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
     * @return a FK4 coordinate reference frame without Eterms
     */    
    private CoordinateReferenceFrame createFK4NOE() {
        final CoordinateReferenceFrame refSystem = new FK4_NO_E();
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
     * <p>
     * The coordinate reference system is found according to the CTYPE1 keyword.
     *
     * @return the coordinate reference system
     * @throws JWcsError The coordinate reference system is not supported
     */
    public Crs getCrs() {
        final Crs crs;
        final CoordinateReferenceFrame refSystem = getReferenceFrame();
        if (hasKeyword("CTYPE1")) {
            String ctype1 = getValueAsString("CTYPE1");
            ctype1 = ctype1.substring(0, ctype1.indexOf('-'));
            switch (ctype1) {
                case "RA":
                    crs = new Equatorial();
                    if (refSystem != null) {
                        ((Equatorial) crs).setCoordinateReferenceFrame(refSystem);
                    }
                    break;
                case "DEC":
                    crs = new Equatorial();
                    if (refSystem != null) {
                        ((Equatorial) crs).setCoordinateReferenceFrame(refSystem);
                    }
                    break;
                case "GLON":
                    crs = new Galactic();
                    break;
                case "GLAT":
                    crs = new Galactic();
                    break;
                case "ELON":
                    crs = new Ecliptic();
                    if (refSystem != null) {
                        ((Ecliptic) crs).setCoordinateReferenceFrame(refSystem);
                    }
                    break;
                case "ELAT":
                    crs = new Ecliptic();
                    if (refSystem != null) {
                        ((Ecliptic) crs).setCoordinateReferenceFrame(refSystem);
                    }
                    break;
                default:
                    throw new JWcsError("The coordinate reference system " + ctype1 + " is not supported");
            }
        } else {
            throw new JWcsError("Cannot find the coordinate reference system.");
        }

        return crs;

    }

    /**
     * Make the initialization of the WCS. By default, you must call
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
     * <p>
     * The supported projection families are the following:
     * <ul>
     * <li>{@link io.github.malapert.jwcs.proj.CylindricalProjection}</li>
     * <li>{@link io.github.malapert.jwcs.proj.ConicProjection}</li>
     * <li>{@link io.github.malapert.jwcs.proj.PolyConicProjection}</li>
     * <li>{@link io.github.malapert.jwcs.proj.ZenithalProjection}</li>
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
     * The computation is realized as follows:
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
        final double cos0 = Math.cos(Math.toRadians(crota));
        final double sin0 = Math.sin(Math.toRadians(crota));
        final double cd11 = cdelt[0] * cos0;
        final double cd12 = Math.abs(cdelt[1]) * Math.signum(cdelt[0]) * sin0;
        final double cd21 = -Math.abs(cdelt[0]) * Math.signum(cdelt[1]) * sin0;
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
     * <p>
     * The CD matrix is created by reading CD matrix or by computing the CD
     * matrix from the CDELT and CROTA.
     * </p>
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
            final double[][] cdTmp = JWcs.computeCdFromCdelt(cdelt, getValueAsDouble(CROTA2));
            result = cdTmp[i - 1][j - 1];
        } else if (hasKeyword(PC11)) {
            final double[][] pc = new double[][]{{getValueAsDouble(PC11), getValueAsDouble(PC12)},
            {getValueAsDouble(PC21), getValueAsDouble(PC22)}};
            final double[] cdelt = new double[]{getValueAsDouble(CDELT1), getValueAsDouble(CDELT2)};
            final double[][] cdTmp = JWcs.pc2cd(pc, cdelt);
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
     * Scale factor. This is used to convert the CRVAL into degree.
     *
     * @param cunit The cunit axis
     * @return the scale factor to apply at CRVAL
     */
    private double convertToDegree(final String cunit) {
        final double cx;
        if (hasKeyword(cunit)) {
            final String unit_lc = cunit.toLowerCase();
            switch (unit_lc) {
                case "acmin":
                    cx = 1 / 60;
                    break;
                case "arcsec":
                    cx = 1 / 3600;
                    break;
                case "mas":
                    cx = 1 / 3600000;
                    break;
                case "rad":
                    cx = 180 / Math.PI;
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
     * <p>
     * Raises a JWcsError when no projection code is found.
     *
     * @return the projection
     * @throws BadProjectionParameterException when the projection parameter is
     * wrong
     */
    protected final Projection createProjection() throws BadProjectionParameterException {
        final String ctype1 = ctype(1);
        final String codeProjection = ctype1.substring(ctype1.lastIndexOf('-') + 1, ctype1.length());
        final double cx = convertToDegree(cunit(1));
        final double cy = convertToDegree(cunit(2));
        final Projection projection = createProjectionFactory(codeProjection, cx, cy);
        setNativeLongitudeOfFiducialPoint(projection);
        setNativeLatitudeOfFiducialPoint(projection);
        setNativeLongitudeOfCelestialPole(projection);
        setNativeLatitudeOfCelestialPole(projection);
        return projection;
    }

    /**
     * Creates a projection based on its projection code
     *
     * @param codeProjection projection code
     * @param cx scale factor along X
     * @param cy scale factor along Y
     * @return the projection
     * @throws BadProjectionParameterException When a bad parameter is provided
     * to the projection
     * @throws JWcsError projection code is not supported
     */
    private Projection createProjectionFactory(final String codeProjection, final double cx, final double cy) throws BadProjectionParameterException {
        final Projection projection;
        switch (codeProjection) {
            case "AIT":
                LOG.log(Level.INFO, "Creates a AIT projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
                projection = new AIT(crval(1) * cx, crval(2) * cy);
                break;
            case "ARC":
                LOG.log(Level.INFO, "Creates a ARC projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
                projection = new ARC(crval(1) * cx, crval(2) * cy);
                break;
            case "AZP":
                projection = createAZPProjection(cx, cy);
                break;
            case "BON":
                LOG.log(Level.INFO, "Creates a AIT projection with (crval1,crval2)=({0},{1}) theta1={2}", new Object[]{crval(1) * cx, crval(2) * cx, getValueAsDouble(PV21, 0)});
                projection = new BON(crval(1) * cx, crval(2) * cy, getValueAsDouble(PV21, 0));
                break;
            case "CAR":
                LOG.log(Level.INFO, "Creates a CAR projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
                projection = new CAR(crval(1) * cx, crval(2) * cy);
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
            case "CYP":
                projection = createCYPProjection(cx, cy);
                break;
            case "MER":
                LOG.log(Level.INFO, "Creates a MER projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
                projection = new MER(crval(1) * cx, crval(2) * cy);
                break;
            case "MOL":
                LOG.log(Level.INFO, "Creates a MOL projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
                projection = new MOL(crval(1) * cx, crval(2) * cy);
                break;
            case "PAR":
                LOG.log(Level.INFO, "Creates a PAR projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
                projection = new PAR(crval(1) * cx, crval(2) * cy);
                break;
            case "PCO":
                LOG.log(Level.INFO, "Creates a PCO projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
                projection = new PCO(crval(1) * cx, crval(2) * cy);
                break;
            case "SFL":
                LOG.log(Level.INFO, "Creates a SFL projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
                projection = new SFL(crval(1) * cx, crval(2) * cy);
                break;
            case "SIN":
                projection = createSINProjection(cy, cy);
                break;
            case "STG":
                LOG.log(Level.INFO, "Creates a STG projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
                projection = new STG(crval(1) * cx, crval(2) * cy);
                break;
            case "SZP":
                projection = createSZPProjection(cx, cy);
                break;
            case "TAN":
                LOG.log(Level.INFO, "Creates a TAN projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
                projection = new TAN(crval(1) * cx, crval(2) * cy);
                break;
            case "ZEA":
                LOG.log(Level.INFO, "Creates a ZEA projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
                projection = new ZEA(crval(1) * cx, crval(2) * cy);
                break;
            case "ZPN":
                projection = createZPNProjection(cx, cy);
                break;
            default:
                throw new JWcsError("code projection : " + codeProjection + " is not supported");
        }
        return projection;
    }

    /**
     * Creates a AZP projection.
     *
     * @param cx scale factor along X
     * @param cy scale factor along Y
     * @return the AZP projection
     * @throws BadProjectionParameterException when a bad parameter is provided
     * to the projection
     */
    private Projection createAZPProjection(final double cx, final double cy) throws BadProjectionParameterException {
        final Projection projection;
        if (hasKeyword(PV21) && hasKeyword(PV22)) {
            final double mu = getValueAsDouble(PV21);
            final double gamma = getValueAsDouble(PV22);
            LOG.log(Level.INFO, "Creates a AIT projection with (crval1,crval2)=({0},{1}) (mu,gamma)=({2},{3})", new Object[]{crval(1) * cx, crval(2) * cx, mu, gamma});
            projection = new AZP(crval(1) * cx, crval(2) * cy, mu, gamma);
        } else {
            LOG.log(Level.INFO, "Creates a AZP projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
            projection = new AZP(crval(1) * cx, crval(2) * cy);
        }
        return projection;
    }

    /**
     * Creates a CYP projection.
     *
     * @param cx scale factor along X
     * @param cy scale factor along Y
     * @return the CYP projection
     * @throws BadProjectionParameterException when a bad parameter is provided
     * to the projection
     */
    private Projection createCYPProjection(final double cx, final double cy) throws BadProjectionParameterException {
        final Projection projection;
        if (hasKeyword(PV21) && hasKeyword(PV22)) {
            LOG.log(Level.INFO, "Creates a CYP projection with (crval1,crval2)=({0},{1}) (mu,lambda)=({2},{3})", new Object[]{crval(1) * cx, crval(2) * cx, getValueAsDouble(PV21, 0), getValueAsDouble(PV22, 0)});
            projection = new CYP(crval(1) * cx, crval(2) * cy, getValueAsDouble(PV21), getValueAsDouble(PV22));
        } else {
            LOG.log(Level.INFO, "Creates a CYP projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
            projection = new CYP(crval(1) * cx, crval(2) * cy);
        }
        return projection;
    }

    /**
     * Creates a SIN projection.
     *
     * @param cx scale factor along X
     * @param cy scale factor along Y
     * @return the SIN projection
     * @throws BadProjectionParameterException when a bad parameter is provided
     * to the projection
     */
    private Projection createSINProjection(final double cx, final double cy) throws BadProjectionParameterException {
        final Projection projection;
        if (hasKeyword(PV21) && hasKeyword(PV22)) {
            LOG.log(Level.INFO, "Creates a SIN projection with (crval1,crval2)=({0},{1}) (ksi,eta)=({2},{3})", new Object[]{crval(1) * cx, crval(2) * cx, getValueAsDouble(PV21, 0), getValueAsDouble(PV22, 0)});
            projection = new SIN(crval(1) * cx, crval(2) * cy, getValueAsDouble(PV21), getValueAsDouble(PV22));
        } else {
            LOG.log(Level.INFO, "Creates a SIN projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
            projection = new SIN(crval(1) * cx, crval(2) * cy);
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
    private Projection createSZPProjection(final double cx, final double cy) throws BadProjectionParameterException {
        final Projection projection;
        if (hasKeyword(PV21) && hasKeyword(PV22) && hasKeyword(PV23)) {
            LOG.log(Level.INFO, "Creates a SZP projection with (crval1,crval2)=({0},{1}) (mu,phic,thetac)=({2},{3},{4})", new Object[]{crval(1) * cx, crval(2) * cx, getValueAsDouble(PV21), getValueAsDouble(PV22), getValueAsDouble(PV23)});
            projection = new SZP(crval(1) * cx, crval(2) * cy, getValueAsDouble(PV21), getValueAsDouble(PV22), getValueAsDouble(PV23));
        } else {
            LOG.log(Level.INFO, "Creates a SZP projection with (crval1,crval2)=({0},{1})", new Object[]{crval(1) * cx, crval(2) * cx});
            projection = new SZP(crval(1) * cx, crval(2) * cy);
        }
        return projection;
    }

    /**
     * Creates ZPN projection
     *
     * @param cx the scale factor along X
     * @param cy the scale factor along Y
     * @return the ZPN projection
     * @throws BadProjectionParameterException when a bad parameter is provided
     * to the projection
     */
    private Projection createZPNProjection(final double cx, final double cy) throws BadProjectionParameterException {
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
     * Sets the Native longitude of the fiducial point to the projection
     *
     * @param projection the projection
     */
    private void setNativeLongitudeOfFiducialPoint(final Projection projection) {
        if (hasKeyword(PV11)) {
            LOG.log(Level.INFO, "Sets phi0 to {0}", getValueAsDouble(PV11));
            projection.setPhi0(getValueAsDouble(PV11));
        }
    }

    /**
     * Sets the Native latitude of the celestial pole to the projection
     *
     * @param projection the projection
     */
    private void setNativeLatitudeOfFiducialPoint(final Projection projection) {
        if (hasKeyword(PV12)) {
            LOG.log(Level.INFO, "Sets theta0 to {0}", getValueAsDouble(PV12));
            projection.setTheta0(getValueAsDouble(PV12));
        }
    }

    /**
     * Sets the Native longitude of the celestial pole to the projection
     *
     * @param projection the projection
     */
    private void setNativeLongitudeOfCelestialPole(final Projection projection) {
        if (!Double.isNaN(lonpole())) {
            LOG.log(Level.INFO, "Sets phip to {0}", lonpole());
            projection.setPhip(Math.toRadians(lonpole()));
        }
    }

    /**
     * Sets the Native latitude of the celestial pole to the projection
     *
     * @param projection the projection
     */
    private void setNativeLatitudeOfCelestialPole(final Projection projection) {
        if (!Double.isNaN(latpole())) {
            LOG.log(Level.INFO, "Sets thetap to {0}", latpole());
            projection.setThetap(Math.toRadians(latpole()));
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
        return this.getProj().inside(Math.toRadians(lon), Math.toRadians(lat));
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
        final double[] coordVal = this.getProj().wcs2projectionPlane(Math.toRadians(longitude), Math.toRadians(latitude));
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
    protected final Projection getProj() {
        return proj;
    }

    /**
     * Sets the projection.
     *
     * @param proj the projection to set
     */
    protected final void setProj(final Projection proj) {
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
     * Sets the CD matrix inverse
     *
     * @param cdInverse the CD matrix inverse to set
     */
    protected final void setCdInverse(final RealMatrix cdInverse) {
        this.cdInverse = cdInverse;
    }
}
