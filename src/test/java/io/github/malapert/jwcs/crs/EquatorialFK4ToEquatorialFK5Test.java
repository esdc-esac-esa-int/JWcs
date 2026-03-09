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

import io.github.malapert.jwcs.datum.CoordinateReferenceFrame;
import io.github.malapert.jwcs.datum.FK4;
import io.github.malapert.jwcs.datum.FK5;

/**
 *
 * @author malapert
 */
public class EquatorialFK4ToEquatorialFK5Test extends ConverterTest {

    /**
     *
     */
    public EquatorialFK4ToEquatorialFK5Test() {
        super();
        final CoordinateReferenceFrame fk5 = new FK5();
        final CoordinateReferenceFrame fk4 = new FK4();
        this.source = new Equatorial(fk4);
        this.target = new Equatorial(fk5);

    }
}
