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
public class EquatorialFK5ToEquatorialFK51950Test extends ConverterTest {
    
    public EquatorialFK5ToEquatorialFK51950Test() throws JWcsException {
        super();
        ReferenceSystemInterface fk5 = new FK5();
        ReferenceSystemInterface fk51950 = new FK5(1950);
        this.source = new Equatorial(fk5);
        this.target = new Equatorial(fk51950);        
        
    }
}