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

import io.github.malapert.jwcs.AbstractJWcs;
import io.github.malapert.jwcs.proj.exception.BadProjectionParameterException;
import io.github.malapert.jwcs.proj.exception.MathematicalSolutionException;
import io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException;
import io.github.malapert.jwcs.utility.AiryFunction;
import io.github.malapert.jwcs.utility.NumericalUtility;
import org.apache.commons.math3.util.FastMath;

/**
 * The Airy projection minimizes the error for the region within \u03B8b
 * (Evenden 1991).
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @see <a href="http://www.atnf.csiro.au/people/mcalabre/WCS/ccs.pdf">
 * "Representations of celestial coordinates in FITS", M. R. Calabretta and E.
 * W. Greisen - page 14</a>
 */
public final class AIR extends AbstractZenithalProjection {

    /**
     * Projection's name.
     */
    private static final String NAME_PROJECTION = "Airy projection";

    /**
     * Projection's description.
     */
    private static final String DESCRIPTION = "\u03B8b=%s";

    /**
     * Default value for \u03B8b sets to 90&deg;.
     */
    public static final double DEFAULT_VALUE_THETHAB = 90;

    /**
     * \u03B8<SUB>b</SUB> value.
     */
    private double thetab;

    /**
     * The Airy function R<sub>\u03B8</sub> to solve in an iterative way when a
     * point is projected.
     *
     * <p>R<sub>\u03B8</sub> = -2 *
     * (ln(cos\u03B6)/tan\u03B6+ln(cos\u03B6<sub>b</sub>)/tan<sup>2</sup>\u03B6<sub>b</sub>*tan\u03B6)
     * with:
     * <ul>
     * <li>\u03B6 = 0.5 * (HALF_PI - \u03B8)</li>
     * <li>\u03B6<sub>b</sub> = 0.5 * (HALF_PI - \u03B8<sub>b</sub>)</li>
     * </ul>
     * @see io.github.malapert.jwcs.utility.AiryFunction
     */
    private final AiryFunction airyFunction;
    
    /**
     * Creates a new AIR projection based on the default \u03B8<sub>b</sub> 
     * and celestial longitude and latitude of the fiducial 
     * point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     *
     * <p>\u03B8<sub>b</sub> is set to {@link AIR#DEFAULT_VALUE_THETHAB}.
     *   
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException \u03B8<sub>b</sub> outside the range ]-90,90]   
     */    
    public AIR() throws BadProjectionParameterException {
        this(FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_THETA0));
    }

    /**
     * Creates a new AIR projection based on the celestial longitude and
     * latitude of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     *
     * <p>\u03B8 is set to {@link AIR#DEFAULT_VALUE_THETHAB}.
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point     
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException \u03B8<sub>b</sub> outside the range ]-90,90]   
     */
    public AIR(final double crval1, final double crval2) throws BadProjectionParameterException {
        this(crval1, crval2, DEFAULT_VALUE_THETHAB);
    }

    /**
     * Creates a new AIR projection based on the celestial longitude and
     * latitude of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     *
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @param thetab \u03B4<sub>b</sub>, projection parameter
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException \u03B8<sub>b</sub> outside the range ]-90,90]
     */
    public AIR(final double crval1, final double crval2, final double thetab) throws BadProjectionParameterException {
        super(crval1, crval2);
        setThetab(FastMath.toRadians(thetab));
        this.airyFunction = new AiryFunction(this.thetab);
    }
    
    /**
     * Checks the validity range of \u03B8<sub>b</sub>.
     * @param thetab \u03B8<sub>b</sub>, projection parameter
     * @throws BadProjectionParameterException \u03B8<sub>b</sub> outside the range ]-90,90]
     */
    private void checkParameter(final double thetab) throws BadProjectionParameterException {
        if (!NumericalUtility.isInInterval(thetab, -NumericalUtility.HALF_PI, false, NumericalUtility.HALF_PI, true)) {
            throw new BadProjectionParameterException(this, "\u03B4<sub>b</sub> is outside the valid range ]-90,90]");
        }
    }

    /**
     * Computes the native spherical coordinates (\u03D5, \u03B8) from the projection plane
     * coordinates (x, y).
     * 
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes the radius : {@link AbstractZenithalProjection#computeRadius(double, double) }</li>
     * <li>solves by an iterative way \u03B8 : {@link AiryFunction}</li>
     * <li>computes \u03D5 : {@link AbstractZenithalProjection#computePhi(double, double, double) }</li>
     * </ul>
     *
     * @param x projection plane coordinate along X
     * @param y projection plane coordinate along Y
     * @return the native spherical coordinates (\u03D5, \u03B8) in radians
     * @throws io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException 
     * No valid solution for (x,y) when an error happens during {@link NumericalUtility#computeFunctionSolution(int, org.apache.commons.math3.analysis.UnivariateFunction, double, double) }
     */    
    @Override
    protected double[] project(final double x, final double y) throws PixelBeyondProjectionException {  
        final double xr = FastMath.toRadians(x);
        final double yr = FastMath.toRadians(y);
        final double radius = this.computeRadius(xr, yr);
        getAiryFunction().setRadius(radius);
        final double theta;
        try {
            theta = NumericalUtility.computeFunctionSolution(1000, getAiryFunction(), -NumericalUtility.HALF_PI, NumericalUtility.HALF_PI);
        } catch(MathematicalSolutionException ex) {
            throw new PixelBeyondProjectionException(this, x, y, ex.getMessage(), true);
        }
        final double phi = computePhi(xr, yr, radius);
        return new double[]{phi, theta};
    }

    /**
     * Computes the projection plane coordinates (x, y) from the native spherical
     * coordinates (\u03D5, \u03B8).
     * 
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes c : {@link AIR#computeC() }</li>
     * <li>computes term1 = ln(cos\u03B6) / tan\u03B6</li>
     * <li>computes the radius : -2 * (term1 + c * tan(\u03B6))</li>
     * <li>computes x : {@link AbstractZenithalProjection#computeX(double, double) }</li>
     * <li>computes y : {@link AbstractZenithalProjection#computeY(double, double) }</li>
     * </ul>     
     *
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param theta the native spherical coordinate (\u03B6) in radians along latitude
     * @return the projection plane coordinates
     * @throws io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException No valid solution for (x,y) when \u03B6=0
     */    
    @Override
    protected double[] projectInverse(final double phi, final double theta) throws PixelBeyondProjectionException  {
        final double c = computeC();
        final double zeta = 0.5 * (NumericalUtility.HALF_PI - theta);
        if (NumericalUtility.equal(zeta, 0)) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(phi), FastMath.toDegrees(theta), false);
        }
        final double term1 = FastMath.log(FastMath.cos(zeta)) / FastMath.tan(zeta);
        final double radius = -2 * (term1 + c * FastMath.tan(zeta));
        final double x = computeX(radius, phi);
        final double y = computeY(radius, phi);
        return new double[]{FastMath.toDegrees(x), FastMath.toDegrees(y)};
    }
    
    /**
     * Computes ln(cos\u03B6<sub>b</sub>) / tan<sup>2</sup>\u03B6<sub>b</sub>.
     * 
     * <p>When cos\u03B6<sub>b</sub> = 0, then 0 is returned<br>
     * When cos\u03B6<sub>b</sub> = 1 or tan\u03B6<sub>b</sub> = 0, then c 
     * approaches to its asymptotic value of -0.5<br>
     * Otherwise computes ln(cos\u03B6<sub>b</sub>) / tan<sup>2</sup>\u03B6<sub>b</sub>
     * 
     * @return c, the result of the computation
     */
    private double computeC() {
        final double zetab = 0.5 * (NumericalUtility.HALF_PI - getThetab());
        final double cos_zetab = FastMath.cos(zetab); 
        final double c;
        if (NumericalUtility.equal(cos_zetab, 0)) {
            c = 0d;
        } else if (NumericalUtility.equal(cos_zetab, 1) || NumericalUtility.equal(FastMath.tan(zetab), 0)) {
            c = -0.5d;
        } else {
            c = FastMath.log(FastMath.cos(zetab))/FastMath.pow(FastMath.tan(zetab), 2);
        }         
        return c;
    }    
    
    @Override
    public boolean inside(final double lon, final double lat) {
        final double raFixed = NumericalUtility.normalizeLongitude(lon);
        final double[] nativeSpherical = computeNativeSpherical(raFixed, lat);
        nativeSpherical[0] = phiRange(nativeSpherical[0]);
        return !NumericalUtility.equal(nativeSpherical[1], -NumericalUtility.HALF_PI);
    }    

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return String.format(DESCRIPTION, NumericalUtility.round(FastMath.toDegrees(this.getThetab())));
    }

    @Override
    public ProjectionParameter[] getProjectionParameters() {
        final ProjectionParameter p1 = new ProjectionParameter("\u03B8b", AbstractJWcs.PV21, new double[]{-90, 90}, DEFAULT_VALUE_THETHAB);
        return new ProjectionParameter[]{p1};
    }

    /**
     * Returns the Airy function to solve.
     *
     * @return the airyFunction     
     */
    private AiryFunction getAiryFunction() {
        return airyFunction;
    }

    /**
     * Returns \u03B8<SUB>b</SUB> in radians.
     *
     * @return the thetab
     */
    public double getThetab() {
        return thetab;
    }

    /**
     * Sets \u03B8<SUB>b</SUB> in radians.
     *
     * @param thetab the thetab to set
     * @throws io.github.malapert.jwcs.proj.exception.BadProjectionParameterException 
     * \u03B8<sub>b</sub> outside the range ]-90,90]
     */
    public void setThetab(final double thetab) throws BadProjectionParameterException {
        checkParameter(thetab);
        this.thetab = thetab;
    }

}
