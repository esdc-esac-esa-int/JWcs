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

import io.github.malapert.jwcs.crs.Ecliptic;
import io.github.malapert.jwcs.crs.Equatorial;
import io.github.malapert.jwcs.datum.FK5;
import io.github.malapert.jwcs.datum.CoordinateReferenceFrame;
import io.github.malapert.jwcs.proj.exception.JWcsException;

/**
 *
 * @author malapert
 */
public class EquatorialToEclipticFK5Test extends ConverterTest {
    
    /**
     *
     * @throws JWcsException
     */
    public EquatorialToEclipticFK5Test() throws JWcsException {
        super();
        final CoordinateReferenceFrame fk5 = new FK5();
        this.source = new Equatorial();
        this.target = new Ecliptic(fk5);        
        
    }
}
