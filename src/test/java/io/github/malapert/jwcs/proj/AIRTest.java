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

import io.github.malapert.jwcs.JWcsFits;
import io.github.malapert.jwcs.proj.exception.JWcsException;
import nom.tam.fits.Fits;
import nom.tam.fits.FitsException;

import java.io.IOException;
import java.util.Objects;

/**
 *
 * @author malapert
 */
public class AIRTest extends AbstractProjectionTest {

    public AIRTest() throws JWcsException, FitsException, IOException {
        super(new JWcsFits(new Fits(
                Objects.requireNonNull(AIRTest.class.getClassLoader().getResource("1904-66_AIR.fits")).toString())));
    }

}
