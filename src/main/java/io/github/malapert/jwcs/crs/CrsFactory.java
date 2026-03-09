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

import io.github.malapert.jwcs.proj.exception.JWcsError;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author malapert
 */
public abstract class CrsFactory {
    /**
     * Logger.
     */
    private static final Logger LOG = Logger.getLogger(CrsFactory.class.getName());
    
    /**
     * Creates a coordinate reference system based on the coordinate system and 
     * a default coordinate reference frame.
     * 
     * <p>The coordinate reference system is built based on ICRS reference frame
     * when the coordinate system is equatorial or ecliptic
     *
     * @param coordinateSystem the coordinate system
     * @return the coordinate reference system
     * @exception JWcsError the coordinate system is not supported
     * @see Ecliptic
     * @see Equatorial
     * @see Galactic
     * @see SuperGalactic
     */    
    public static AbstractCrs create(final AbstractCrs.CoordinateReferenceSystem coordinateSystem) {
        AbstractCrs crs;
        LOG.log(Level.INFO, "Get sky system {0}", new Object[]{coordinateSystem.name()});
        switch (coordinateSystem) {
            case ECLIPTIC:
                crs = new Ecliptic();
                break;
            case EQUATORIAL:
                crs = new Equatorial();
                break;
            case GALACTIC:
                crs = new Galactic();
                break;
            case SUPER_GALACTIC:
                crs = new SuperGalactic();
                break;
            default:
                throw new JWcsError(coordinateSystem + " not supported as coordinate reference system");
        }
        return crs;        
    }
    
}
