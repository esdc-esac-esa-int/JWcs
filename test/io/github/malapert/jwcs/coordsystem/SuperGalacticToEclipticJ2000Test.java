/*
 * Copyright (C) 2016 malapert
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package io.github.malapert.jwcs.coordsystem;

import io.github.malapert.jwcs.proj.exception.JWcsException;

/**
 *
 * @author malapert
 */
public class SuperGalacticToEclipticJ2000Test extends ConverterTest {
    
    public SuperGalacticToEclipticJ2000Test() throws JWcsException {
        super();
        ReferenceSystemInterface j2000 = new J2000();
        this.source = new SuperGalactic();
        this.target = new Ecliptic(j2000);        
        
    }
}