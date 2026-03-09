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
package io.github.malapert.jwcs.utility;

import io.github.malapert.jwcs.proj.exception.MathematicalSolutionException;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.util.FastMath;

/**   
 * The Airy function R<sub>\u03B8</sub> to solve in an iterative way when a
 * point is projected.
 *
 * <p>R<sub>\u03B8</sub> + 2 *
 * ( ln(cos\u03B6) / tan\u03B6 + ln(cos\u03B6<sub>b</sub>) / tan<sup>2</sup>\u03B6<sub>b</sub> * tan\u03B6) = 0
 * with:
 * <ul>
 * <li>\u03B6 = 0.5 * (HALF_PI - \u03B8)</li>
 * <li>\u03B6<sub>b</sub> = 0.5 * (HALF_PI - \u03B8<sub>b</sub>)</li>
 * </ul>     
 * 
 * <p>Special cases must be handled for ln(cos\u03B6<sub>b</sub>) / tan<sup>2</sup>\u03B6<sub>b</sub>,
 * this is managed to {@link AiryFunction#computeC() }
 * 
 * <p>When \u03B8<sub>b</sub>
 * 
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 */
public class AiryFunction implements UnivariateFunction {
    
    /**
     * \u03B8<sub>b</sub>.
     */
    private final double thetab;
    
    /**
     * Radius.
     */
    private double radius;
    
    /**
     * Initialize the function with \u03B8<sub>b</sub> with radius = 0.
     * 
     * <p>\u03B8<sub>b</sub> is a fixed parameter. Only the radius changes
     * according to (x,y) values.
     * 
     * @param thetab \u03B8<sub>b</sub>.
     */
    public AiryFunction(final double thetab) {
        this.thetab = thetab;
        this.radius = 0;
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
        if (NumericalUtility.equalValues(cos_zetab, 0)) {
            c = 0d;
        } else if (NumericalUtility.equalValues(cos_zetab, 1) || NumericalUtility.equalValues(FastMath.tan(zetab), 0)) {
            c = -0.5d;
        } else {
            c = FastMath.log(FastMath.cos(zetab))/FastMath.pow(FastMath.tan(zetab), 2);
        }         
        return c;
    }

    @Override
    public double value(final double theta) {
        final double c = computeC();
        final double zeta = 0.5 * (NumericalUtility.HALF_PI - theta);        
        if (NumericalUtility.equalValues(zeta, 0)) {
            throw new MathematicalSolutionException("zeta cannot be 0");
        }
        final double lncZeta = FastMath.log(FastMath.cos(zeta));
        final double tanZeta = FastMath.tan(zeta);
        return 2 * (lncZeta / tanZeta + c * tanZeta) + getRadius();       
    }

    /**
     * Returns \u03B8<sub>b</sub>.
     * @return the thetab
     */
    public double getThetab() {
        return thetab;
    }

    /**
     * Returns the radius.
     * @return the radius
     */
    public double getRadius() {
        return radius;
    }

    /**
     * Sets the radius.
     * @param radius the radius to set
     */
    public void setRadius(final double radius) {
        this.radius = radius;
    }
    
}
