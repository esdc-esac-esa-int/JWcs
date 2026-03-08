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

import io.github.malapert.jwcs.datum.J2000;

/**
 *
 * @author malapert
 */
public class EquatorialJ20001950EquatorialJ20002001Test extends ConverterTest {

    public EquatorialJ20001950EquatorialJ20002001Test() {
        J2000 j1 = new J2000();
        j1.setEquinox(1950);
        J2000 j2 = new J2000();
        j2.setEquinox(2001);
        this.source = new Equatorial(j1);
        this.target = new Equatorial(j2);
    }

}
