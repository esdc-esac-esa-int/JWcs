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

/**
 * Projection Exception.
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 */
public class ProjectionException extends JWcsException {
    private final static long serialVersionUID = 3428780035222139573L;
    
    /**
     * Class where the exception happens.
     */
    private final AbstractProjection projectionName;
    
    public ProjectionException(final AbstractProjection projectionName) {
        super();  
        this.projectionName = projectionName;
    }

    /**
     * Creates a ProjectionException based on the projection class and a message.
     * @param projectionName projection class
     * @param message message
     */
    public ProjectionException(final AbstractProjection projectionName, final String message) {
        super(message);
        this.projectionName = projectionName;        
    }
    
    /**
     * Returns the projection name.
     * @return the projection name
     */
    public final AbstractProjection getProjection() {
        return this.projectionName;
    }

    @Override
    public String toString() {
        return this.projectionName.getName()+" - "+this.getMessage();
    }           
}
