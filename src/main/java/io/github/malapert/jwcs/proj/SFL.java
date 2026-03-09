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

import io.github.malapert.jwcs.utility.NumericalUtility;
import java.util.logging.Level;
import org.apache.commons.math3.util.FastMath;

/**
 * Sanson-Flamsteed.
 * 
 * <p>Bonne's projection reduces to the pseudocylindrical Sanson-Flamsteed16 
 * projection when theta1 = 0. Parallels are equispaced and projected at 
 * their true length which makes it an equal area projection.
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class SFL extends AbstractCylindricalProjection {
    
    /**
     * Projection's name.
     */
    private static final String NAME_PROJECTION = "Sanson-Flamsteed";
    
    /**
     * Projection's description.
     */
    private static final String DESCRIPTION = "no limits";

   /**
     * Constructs a SFL projection based on the default celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     */    
    public SFL() {
        this(FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_PHI0), FastMath.toDegrees(AbstractCylindricalProjection.DEFAULT_THETA0));
    }
    
   /**
     * Constructs a SFL projection based on the celestial longitude and latitude
     * of the fiducial point (\u03B1<sub>0</sub>, \u03B4<sub>0</sub>).
     * 
     * @param crval1 Celestial longitude \u03B1<sub>0</sub> in degrees of the
     * fiducial point
     * @param crval2 Celestial longitude \u03B4<sub>0</sub> in degrees of the
     * fiducial point
     */
    public SFL(final double crval1, final double crval2) {
        super(crval1, crval2);
        LOG.log(Level.FINER, "INPUTS[Deg] (crval1,crval2)=({0},{1})", new Object[]{crval1,crval2});                                        
    }

    @Override
    protected double[] project(final double x, final double y) {        
        final double theta = FastMath.toRadians(y);
        final double cosTheta = FastMath.cos(theta);
        final double phi;
        if(NumericalUtility.equalValues(cosTheta, 0)) {
            phi = 0;
        } else {
            phi = FastMath.toRadians(x) / FastMath.cos(theta);
        }
        return new double[] { phi, theta};
    }

    @Override
    protected double[] projectInverse(final double phi, final double theta) {
        final double x = FastMath.toDegrees(phi * FastMath.cos(theta));
        final double y = FastMath.toDegrees(theta);
        return new double[] { x, y};
    }

    @Override
    public String getName() {
        return NAME_PROJECTION;
    }

    @Override
    public String getDescription() {
        return DESCRIPTION;
    }
}
