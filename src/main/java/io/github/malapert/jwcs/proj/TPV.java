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

import io.github.malapert.jwcs.proj.exception.BadProjectionParameterException;
import io.github.malapert.jwcs.proj.exception.JWcsError;
import io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException;
import io.github.malapert.jwcs.proj.exception.ProjectionException;
import io.github.malapert.jwcs.utility.NumericalUtility;

import static io.github.malapert.jwcs.utility.NumericalUtility.HALF_PI;

import java.util.Arrays;
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * Gnomonic.
 *
 * <p>The zenithal perspective projection with mu = 0, the gnomonic projection, is
 * widely used in optical astronomy and was given its own code within the AIPS
 * convention, namely TAN. The implementation in this class adds the distortion 
 * parameters PVi_j.
 *
 * <p>Ref : http://www.atnf.csiro.au/people/mcalabre/WCS/ccs.pdf page 12
 *
 * @author talonso - Tomas.AlonsoAlbi@esa.int
 * @version 1.0
 */
public class TPV extends AbstractZenithalProjection {

    /**
     * Projection's name.
     */
    private static final String NAME_PROJECTION = "Gnomonic with correction terms";
    
    /**
     * Projection's description.
     */
    private static final String DESCRIPTION = "poly=(%s)";     

    /**
     * Projection parameters for longitude.
     */
    private double[] a = new double[40];
    /**
     * Projection parameters for latitude.
     */
    private double[] b = new double[40];
    
   /**
     * Constructs a TPV projection based on the default celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>). Projection 
     * parameters are set to null, so this is equivalent to TAN.
     * @throws BadProjectionParameterException when a parameter projection is wrong
     */    
    public TPV() throws BadProjectionParameterException {
        this(FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractZenithalProjection.DEFAULT_THETA0), null, null);
    }

   /**
     * Constructs a TPV projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     * 
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     * @param a projection parameters for longitude
     * @param b projection parameters for latitude
     * @throws BadProjectionParameterException when a parameter projection is wrong
     */
    public TPV(final double crval1, final double crval2, double[] a, double[] b) throws BadProjectionParameterException {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2)=({0},{1})", new Object[]{crval1,crval2});                                
        setPv(a, b);
        check();
    }

    /**
     * Checks validity of parameters.
     *
     * @throws JWcsError Non-standard PVi_1 and/or PVi_2 values
     * @throws ArrayIndexOutOfBoundsException Need at least 10 projection
     * parameters
     */
    private void check() {
        if (!NumericalUtility.equal(getPhi0(), 0) || !NumericalUtility.equal(getTheta0(), HALF_PI)) {
            throw new JWcsError("Non-standard PVi_1 and/or PVi_2 values");
        }
    }

    /**
     * Computes the native spherical coordinates (\u03D5, \u03B8) from the projection plane
     * coordinates (x, y).
     * 
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes radius : {@link TPV#computeRadius(double, double) }</li>
     * <li>computes \u03D5 : {@link AbstractZenithalProjection#computePhi(double, double, double) }</li>      
     * <li>computes \u03B8 : arg(radius, 1)</li>
     * </ul>
     * 
     * @param x projection plane coordinate along X
     * @param y projection plane coordinate along Y
     * @return the native spherical coordinates (\u03D5, \u03B8) in radians     
     */      
    @Override
    public double[] project(final double x, final double y) {
        final double xr = FastMath.toRadians(x);
        final double yr = FastMath.toRadians(y);
        final double r_theta = computeRadius(xr, yr);
        
        final double[] xy = applyPoly(xr, yr, r_theta);
        final double r_theta_corrected = computeRadius(xy[0], xy[1]);
        
        final double phi = computePhi(xy[0], xy[1], r_theta_corrected);       
        final double theta = NumericalUtility.aatan2(1, r_theta_corrected);        
        return new double[] {phi, theta};
    }

    /**
     * Computes the projection plane coordinates (x, y) from the native spherical
     * coordinates (\u03D5, \u03B8).
     *
     * <p>The algorithm to make this projection is the following:
     * <ul>
     * <li>computes radius : cos\u03B8 / sin\u03B8</li>
     * <li>computes x : {@link AbstractZenithalProjection#computeX(double, double) }</li>
     * <li>computes y : {@link AbstractZenithalProjection#computeY(double, double) }</li>
     * </ul>
     * 
     * @param phi the native spherical coordinate (\u03D5) in radians along longitude
     * @param theta the native spherical coordinate (\u03B8) in radians along latitude
     * @return the projection plane coordinates
     * @throws io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException No valid solution for (\u03D5, \u03B8)
     * @throws io.github.malapert.jwcs.proj.exception.ProjectionException Inexistent solution for (\u03D5, \u03B8), or did
     * not converge
     */     
    @Override
    public double[] projectInverse(final double phi, final double theta) throws PixelBeyondProjectionException, 
        ProjectionException {        
        final double s = FastMath.sin(theta);
        if (NumericalUtility.equal(s, 0)) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(phi), FastMath.toDegrees(theta), false);
        }
        final double r_theta = FastMath.cos(theta) / s;
        final double x = computeX(r_theta, phi);
        final double y = computeY(r_theta, phi);
        
        double[] xy = applyPolyInverse(x, y);
        
        return new double[] {FastMath.toDegrees(xy[0]), FastMath.toDegrees(xy[1])};
    }       

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }
    
    @Override
    public String getDescription() {
        return DESCRIPTION;
    }
    
    @Override
    public boolean inside(final double lon, final double lat) {
        final double raFixed = NumericalUtility.normalizeLongitude(lon);
        final double[] nativeSpherical = computeNativeSpherical(raFixed, lat);
        nativeSpherical[0] = phiRange(nativeSpherical[0]);
        final boolean result = NumericalUtility.equal(nativeSpherical[1], 0);
        return !result && super.inside(lon, lat);
    }        

    @Override
    public ProjectionParameter[] getProjectionParameters() {
        return new ProjectionParameter[]{};
    }

    /**
     * Returns a copy of pv for longitude.
     *
     * @return the pv for longitude
     */
    public double[] getPvLon() {
        return a.clone();
    }

    /**
     * Returns a copy of pv for latitude.
     *
     * @return the pv for latitude
     */
    public double[] getPvLat() {
        return b.clone();
    }

    /**
     * Sets pv.
     *
     * @param a the pv to set for longitude
     * @param b the pv to set for latitude
     * @throws BadProjectionParameterException when a parameter projection is wrong
     */
    public void setPv(final double[] a, final double[] b) throws BadProjectionParameterException {
        String errMsg = "Projection parameters are of different size, or one is null and not the other";
        if (a == null && b != null) 
            throw new BadProjectionParameterException(this, errMsg);
        if (a != null && b == null) 
            throw new BadProjectionParameterException(this, errMsg);
        if (a != null && a.length != b.length) 
            throw new BadProjectionParameterException(this, errMsg);
        
        if (a == null) {
            this.a = new double[0];
        } else {
            this.a = Arrays.copyOf(a, a.length);
        }
        
        if (b == null) {
            this.b = new double[0];
        } else {
            this.b = Arrays.copyOf(b, b.length);
        }
    }
    
    // Implementation below to apply the distortion terms follows the C code at
    // https://github.com/Starlink/ast/blob/master/wcslib/tpn.c
    private double[] applyPoly(double x, double y, double r) {
        if (a == null || a.length == 0) return new double[] {x, y};
        
        double x2 = x*x;
        double xy = x*y;
        double y2 = y*y;

        double r2 = x2 + y2;

        double x3 = x2*x;
        double x2y = x2*y;
        double xy2 = x*y2;
        double y3 = y*y2;

        double r3 = r*r2;

        double x4 = x3*x;
        double x3y = x3*y;
        double x2y2 = x2*y2;
        double xy3 = x*y3;
        double y4 = y*y3;

        double x5 = x4*x;
        double x4y = x4*y;
        double x3y2 = x3*y2;
        double x2y3 = x2*y3;
        double xy4 = x*y4;
        double y5 = y*y4;

        double r5 = r3*r2;

        double x6 = x5*x;
        double x5y = x5*y;
        double x4y2 = x4*y2;
        double x3y3 = x3*y3;
        double x2y4 = x2*y4;
        double xy5 = x*y5;
        double y6 = y*y5;

        double x7 = x6*x;
        double x6y = x6*y;
        double x5y2 = x5*y2;
        double x4y3 = x4*y3;
        double x3y4 = x3*y4;
        double x2y5 = x2*y5;
        double xy6 = x*y6;
        double y7 = y*y6;

        double r7 = r5*r2;

        double xi =   
               a[0]       + a[1]*x     + a[2]*y     + a[3]*r     + a[4]*x2
             + a[5]*xy    + a[6]*y2    + a[7]*x3    + a[8]*x2y   + a[9]*xy2
             + a[10]*y3   + a[11]*r3   + a[12]*x4   + a[13]*x3y  + a[14]*x2y2
             + a[15]*xy3  + a[16]*y4   + a[17]*x5   + a[18]*x4y  + a[19]*x3y2
             + a[20]*x2y3 + a[21]*xy4  + a[22]*y5   + a[23]*r5   + a[24]*x6
             + a[25]*x5y  + a[26]*x4y2 + a[27]*x3y3 + a[28]*x2y4 + a[29]*xy5
             + a[30]*y6   + a[31]*x7   + a[32]*x6y  + a[33]*x5y2 + a[34]*x4y3
             + a[35]*x3y4 + a[36]*x2y5 + a[37]*xy6  + a[38]*y7   + a[39]*r7;

        double eta =  
               b[0]       + b[1]*y     + b[2]*x     + b[3]*r     + b[4]*y2
             + b[5]*xy    + b[6]*x2    + b[7]*y3    + b[8]*xy2   + b[9]*x2y
             + b[10]*x3   + b[11]*r3   + b[12]*y4   + b[13]*xy3  + b[14]*x2y2
             + b[15]*x3y  + b[16]*x4   + b[17]*y5   + b[18]*xy4  + b[19]*x2y3
             + b[20]*x3y2 + b[21]*x4y  + b[22]*x5   + b[23]*r5   + b[24]*y6
             + b[25]*xy5  + b[26]*x2y4 + b[27]*x3y3 + b[28]*x4y2 + b[29]*x5y
             + b[30]*x6   + b[31]*y7   + b[32]*xy6  + b[33]*x2y5 + b[34]*x3y4
             + b[35]*x4y3 + b[36]*x5y2 + b[37]*x6y  + b[38]*x7   + b[39]*r7;
        return new double[] {xi, eta};
    }
    
    private double[] applyPolyInverse(double xi, double eta) throws ProjectionException {
        if (a == null || a.length == 0) return new double[] {xi, eta};
        
        /* Initial guess: linear solution assuming a3,... and b3,... are zero. */
        double denom = a[1]*b[1] - a[2]*b[2];
        double x;
        double y;
        if (denom != 0.0) {
           x = (xi*b[1] - eta*a[2] - a[0]*b[1] + b[0]*a[2]) / denom;
           y = -(xi*b[2] - eta*a[1] - a[0]*b[2] + b[0]*a[1]) / denom;
        } else {
           x = a[0];
           if (a[1] != 0.0)
              x = (xi - a[0]) / a[1];

           y = b[0];
           if (b[1] != 0.0)
              y = (eta - b[0]) / b[1];
        }

        /* Iterate up to niter times, or until the required relative accuracy is
         achieved. Otherwise, find the best values. */
        double tol = 1.0E-5;
        int niter = 200;
        double minD = -1;
        double bestX = x;
        double bestY = y;
        for (int i = 0; i < niter; i++) {
            
            /* Get required products of the current x and y values */
            double x2 = x*x;
            double xy = x*y;
            double y2 = y*y;

            double r2 = x2 + y2;
            double r = FastMath.sqrt(r2);

            double x3 = x2*x;
            double x2y = x2*y;
            double xy2 = x*y2;
            double y3 = y*y2;

            double r3 = r*r2;

            double x4 = x3*x;
            double x3y = x3*y;
            double x2y2 = x2*y2;
            double xy3 = x*y3;
            double y4 = y*y3;

            double x5 = x4*x;
            double x4y = x4*y;
            double x3y2 = x3*y2;
            double x2y3 = x2*y3;
            double xy4 = x*y4;
            double y5 = y*y4;

            double r5 = r3*r2;

            double x6 = x5*x;
            double x5y = x5*y;
            double x4y2 = x4*y2;
            double x3y3 = x3*y3;
            double x2y4 = x2*y4;
            double xy5 = x*y5;
            double y6 = y*y5;

            double x7 = x6*x;
            double x6y = x6*y;
            double x5y2 = x5*y2;
            double x4y3 = x4*y3;
            double x3y4 = x3*y4;
            double x2y5 = x2*y5;
            double xy6 = x*y6;
            double y7 = y*y6;

            double r7 = r5*r2;

            /* Get the xi and eta models corresponding to the current x and y values */
            double f =  
                  a[0]       + a[1]*x     + a[2]*y     + a[3]*r     + a[4]*x2
                + a[5]*xy    + a[6]*y2    + a[7]*x3    + a[8]*x2y   + a[9]*xy2
                + a[10]*y3   + a[11]*r3   + a[12]*x4   + a[13]*x3y  + a[14]*x2y2
                + a[15]*xy3  + a[16]*y4   + a[17]*x5   + a[18]*x4y  + a[19]*x3y2
                + a[20]*x2y3 + a[21]*xy4  + a[22]*y5   + a[23]*r5   + a[24]*x6
                + a[25]*x5y  + a[26]*x4y2 + a[27]*x3y3 + a[28]*x2y4 + a[29]*xy5
                + a[30]*y6   + a[31]*x7   + a[32]*x6y  + a[33]*x5y2 + a[34]*x4y3
                + a[35]*x3y4 + a[36]*x2y5 + a[37]*xy6  + a[38]*y7   + a[39]*r7;

            double g = 
                  b[0]       + b[1]*y     + b[2]*x     + b[3]*r     + b[4]*y2
                + b[5]*xy    + b[6]*x2    + b[7]*y3    + b[8]*xy2   + b[9]*x2y
                + b[10]*x3   + b[11]*r3   + b[12]*y4   + b[13]*xy3  + b[14]*x2y2
                + b[15]*x3y  + b[16]*x4   + b[17]*y5   + b[18]*xy4  + b[19]*x2y3
                + b[20]*x3y2 + b[21]*x4y  + b[22]*x5   + b[23]*r5   + b[24]*y6
                + b[25]*xy5  + b[26]*x2y4 + b[27]*x3y3 + b[28]*x4y2 + b[29]*x5y
                + b[30]*x6   + b[31]*y7   + b[32]*xy6  + b[33]*x2y5 + b[34]*x3y4
                + b[35]*x4y3 + b[36]*x5y2 + b[37]*x6y  + b[38]*x7   + b[39]*r7;

            /* Partial derivative of xi wrt x... */
            double fx = a[1] + a[3]*( (r!=0.0)?(x/r):0.0 ) + 2*a[4]*x +
                a[5]*y       + 3*a[7]*x2    + 2*a[8]*xy    + a[9]*y2 +
                3*a[11]*r*x  + 4*a[12]*x3   + 3*a[13]*x2y  + 2*a[14]*xy2  +
                a[15]*y3     + 5*a[17]*x4   + 4*a[18]*x3y  + 3*a[19]*x2y2 +
                2*a[20]*xy3  + a[21]*y4     + 5*a[23]*r3*x + 6*a[24]*x5  +
                5*a[25]*x4y  + 4*a[26]*x3y2 + 3*a[27]*x2y3 + 2*a[28]*xy4  +
                a[29]*y5     + 7*a[31]*x6   + 6*a[32]*x5y  + 5*a[33]*x4y2 +
                4*a[34]*x3y3 + 3*a[35]*x2y4 + 2*a[36]*xy5  + a[37]*y6 +
                7*a[39]*r5*x;

            /* Partial derivative of xi wrt y... */
            double fy = a[2] + a[3]*( (r!=0.0)?(y/r):0.0 ) + a[5]*x +
                2*a[6]*y     + a[8]*x2      + 2*a[9]*xy    + 3*a[10]*y2 +
                3*a[11]*r*y  + a[13]*x3     + 2*a[14]*x2y  + 3*a[15]*xy2 +
                4*a[16]*y3   + a[18]*x4     + 2*a[19]*x3y  + 3*a[20]*x2y2 +
                4*a[21]*xy3  + 5*a[22]*y4   + 5*a[23]*r3*y + a[25]*x5 +
                2*a[26]*x4y  + 3*a[27]*x3y2 + 4*a[28]*x2y3 + 5*a[29]*xy4 +
                6*a[30]*y5   + a[32]*x6     + 2*a[33]*x5y  + 3*a[34]*x4y2 +
                4*a[35]*x3y3 + 5*a[36]*x2y4 + 6*a[37]*xy5  + 7*a[38]*y6 +
                7*a[39]*r5*y;

            /* Partial derivative of eta wrt x... */
            double gx = b[2] + b[3]*( (r!=0.0)?(x/r):0.0 ) + b[5]*y +
                2*b[6]*x     + b[8]*y2      + 2*b[9]*xy    + 3*b[10]*x2 +
                3*b[11]*r*x  + b[13]*y3     + 2*b[14]*xy2  + 3*b[15]*x2y +
                4*b[16]*x3   + b[18]*y4     + 2*b[19]*xy3  + 3*b[20]*x2y2 +
                4*b[21]*x3y  + 5*b[22]*x4   + 5*b[23]*r3*x + b[25]*y5 +
                2*b[26]*xy4  + 3*b[27]*x2y3 + 4*b[28]*x3y2 + 5*b[29]*x4y +
                6*b[30]*x5   + b[32]*y6     + 2*b[33]*xy5  + 3*b[34]*x2y4 +
                4*b[35]*x3y3 + 5*b[36]*x4y2 + 6*b[37]*x5y  + 7*b[38]*x6 +
                7*b[39]*r5*x;

            /* Partial derivative of eta wrt y... */
            double gy = b[1] + b[3]*( (r!=0.0)?(y/r):0.0 ) + 2*b[4]*y +
                b[5]*x       + 3*b[7]*y2    + 2*b[8]*xy    + b[9]*x2 +
                3*b[11]*r*y  + 4*b[12]*y3   + 3*b[13]*xy2  + 2*b[14]*x2y  +
                b[15]*x3     + 5*b[17]*y4   + 4*b[18]*xy3  + 3*b[19]*x2y2 +
                2*b[20]*x3y  + b[21]*x4     + 5*b[23]*r3*y + 6*b[24]*y5  +
                5*b[25]*xy4  + 4*b[26]*x2y3 + 3*b[27]*x3y2 + 2*b[28]*x4y  +
                b[29]*x5     + 7*b[31]*y6   + 6*b[32]*xy5  + 5*b[33]*x2y4 +
                4*b[34]*x3y3 + 3*b[35]*x4y2 + 2*b[36]*x5y  + b[37]*x6     +
                7*b[39]*r5*y;

            /* Calculate new x and y values. */
            f = f - xi;
            g = g - eta;
            double dx = ( (-f*gy) + (g*fy) ) / ( (fx*gy) - (fy*gx) );
            double dy = ( (-g*fx) + (f*gx) ) / ( (fx*gy) - (fy*gx) );
            x += dx;
            y += dy;
            
            /* Check if convergence has been achieved. */
            if (FastMath.abs(dx) <= tol*FastMath.abs(x) && FastMath.abs(dy) <= tol*FastMath.abs(y) ) 
                return new double[] {x, y};

            // Update best output found so far
            double d = FastMath.hypot(dx, dy);
            if (d < minD || minD < 0) {
                minD = d;
                bestX = x;
                bestY = y;
            }
        }

        if (maximumPixelErrorWhenInvertingDistortionParameters > 0 && FastMath.toDegrees(minD) < maximumPixelErrorWhenInvertingDistortionParameters) {
            LOG.log(Level.FINER, "TPV projection - inversion of projection parameters did not converge - returning best point found");
            return new double[] {bestX, bestY};
        }

        // IMPORTANT:
        // Note that due to the distortion parameters some positions could not exit in the TPV projection, 
        // so the numerical inversion will not converge. It is responsibility of the fits creators to include
        // consistent headers to ensure the projection can be inverted in all points. See Calabretta et al. 2004, 
        // "Representations of distortions in FITS world coordinate systems", page 3
        throw new ProjectionException(this, "TPV projection - inversion of distortion parameters did not converge");
    }
    
    /**
     * Maximum allowed error in pixels when inverting the PV distortion parameters numerically.
     * Default is 0, so in case the inversion does not converge an exception is thrown even if 
     * the process almost converged.
     */
    private static double maximumPixelErrorWhenInvertingDistortionParameters = 0;
    
    /**
     * Sets the maximum pixel error when inverting the projection. Default is 0. This is used 
     * to increase tolerance in case the numerical process does not converge. Set this value to 
     * a fraction of pixel to still return an acceptable output avoiding the exception in some 
     * cases. But note this may be strictly incorrect because the inversion may not exist at 
     * some specific points.
     * @param d
     */
    public static void setMaximumPixelErrorWhenInvertingDistortionParameters(double d) {
        maximumPixelErrorWhenInvertingDistortionParameters = d;
    }
    
    /**
     * Returns the maximum pixel error when inverting the projection. Default is 0.
     * @return Maximum pixel error for the numerical process to invert the projection.
     */
    public static double getMaximumPixelErrorWhenInvertingDistortionParameters() {
        return maximumPixelErrorWhenInvertingDistortionParameters;
    }
}
