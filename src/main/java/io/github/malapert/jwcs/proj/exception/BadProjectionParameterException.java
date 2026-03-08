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

package io.github.malapert.jwcs.proj.exception;

import io.github.malapert.jwcs.proj.AbstractProjection;
import java.util.logging.Level;

/**
 * Bad AbstractProjection Parameter Exception.
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 */
public class BadProjectionParameterException extends ProjectionException{
    private static final long serialVersionUID = -787228345305370079L;
    
    /**
     * Creates a new Exception when the projection parameter are wrong.
     * @param projectionName projection name
     * @param message message
     */
    public BadProjectionParameterException(final AbstractProjection projectionName, final String message) {
        super(projectionName, message);
        projectionName.getLogger().log(Level.SEVERE, "{0}- Bad projection parameter for {1}", new Object[]{getProjection().getClass().getName(), getMessage()});
    }

    @Override
    public String toString() {
        return this.getProjection().getClass().getName()+"- Bad projection parameter for "+getMessage();
    }
    
    
    
}
