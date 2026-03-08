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

import io.github.malapert.jwcs.proj.exception.PixelBeyondProjectionException;
import io.github.malapert.jwcs.utility.NumericalUtility;
import static io.github.malapert.jwcs.utility.NumericalUtility.HALF_PI;
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * Mercator.
 * 
 * <p>Since the meridians and parallels of all cylindrical projections
 * intersect at right angles the requirement for conformality.
 * This projection has been widely used in navigation since it has the 
 * property that lines of constant bearing (known asrhumb lines orloxodromes)
 * are projected as straight lines. 
 * This is a direct result of its conformality and the fact that its meridians 
 * do not converge.
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class MER extends AbstractCylindricalProjection {
    
    /**
     * Projection's name.
     */
    private static final String NAME_PROJECTION = "Mercator";
    
    /**
     * Projection's description.
     */
    private static final String DESCRIPTION = "no limits";

   /**
     * Constructs a MER projection based on the default celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     */    
    public MER() {
        this(FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_THETA0));
    }
    
   /**
     * Constructs a MER projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     * 
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     */
    public MER(final double crval1, final double crval2) {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2)=({0},{1})", new Object[]{crval1,crval2});                
    }

    @Override
    protected double[] project(final double x, final double y) {
        final double xr = FastMath.toRadians(x);
        final double yr = FastMath.toRadians(y);
        final double phi = xr;
        final double theta = 2*FastMath.atan(FastMath.exp(yr)) - HALF_PI;
        return new double[] { phi, theta};
    }

    @Override
    protected double[] projectInverse(final double phi, final double theta) throws PixelBeyondProjectionException  {
        final double x = phi;
        if (NumericalUtility.equal(FastMath.abs(theta), HALF_PI)) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(phi), FastMath.toDegrees(theta), false);            
        }
        final double angle = (HALF_PI + theta) * 0.5d;
        if(NumericalUtility.equal(FastMath.abs(angle), HALF_PI)) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(phi), FastMath.toDegrees(theta), false);      
        }
        final double d = FastMath.tan(angle);
        if (d<0 || NumericalUtility.equal(d, 0)) {
            throw new PixelBeyondProjectionException(this, FastMath.toDegrees(phi), FastMath.toDegrees(theta), false);     
        }
        final double y = FastMath.log(d);
        return new double[] { FastMath.toDegrees(x), FastMath.toDegrees(y)};
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
        return super.inside(lon, lat) && !NumericalUtility.equal(FastMath.abs(lat), HALF_PI);   
    }
}
