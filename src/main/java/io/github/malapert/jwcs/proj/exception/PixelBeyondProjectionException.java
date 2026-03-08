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
import java.text.MessageFormat;
import java.util.logging.Level;

/**
 * Pixel Beyond Projection Exception.
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 */
public class PixelBeyondProjectionException extends ProjectionException {

    private static final long serialVersionUID = -3719985099230583709L;

    /**
     * x coordinate.
     */
    private final double x;
    
    /**
     * y coordinate.
     */
    private final double y;
    
    /**
     * Specifies if the coordinates are planes or spherical.
     */
    private final boolean isPlaneCoordinate;

    /**
     * Creates a PixelBeyondProjectionException based on the pixel coordinate.
     *
     * @param projectionName AbstractProjection
     * @param x x coordinate of the pixel
     * @param y y coordinate of the pixel
     * @param isPlaneCoordinate True when the error is related to plane coordinate
     */
    public PixelBeyondProjectionException(final AbstractProjection projectionName, final double x, final double y, final boolean isPlaneCoordinate) {
        super(projectionName);
        this.x = x;
        this.y = y;
        this.isPlaneCoordinate = isPlaneCoordinate;        
        if (isPlaneCoordinate) {
            getProjection().getLogger().log(Level.FINE, "{0} - Solution not defined for (x,y) = ({1},{2})", new Object[]{this.getProjection().getClass().getName(), x, y});
        } else {
            getProjection().getLogger().log(Level.FINE, "{0} - Solution not defined for (phi,theta) = ({1},{2})", new Object[]{this.getProjection().getClass().getName(), x, y});    
        }
    }

    /**
     * Creates a PixelBeyondProjectionException based on the pixel coordinate
     * and a message.
     *
     * @param projectionName AbstractProjection
     * @param x x coordinate of the pixel
     * @param y y coordinate of the pixel
     * @param message message
     * @param isPlaneCoordinate True when the error is related to plane coordinate
     */
    public PixelBeyondProjectionException(final AbstractProjection projectionName, final double x, final double y, final String message, final boolean isPlaneCoordinate) {
        super(projectionName, message);
        this.x = x;
        this.y = y;
        this.isPlaneCoordinate = isPlaneCoordinate;
        if (isPlaneCoordinate) {
            getProjection().getLogger().log(Level.FINE, "{0} - Solution not defined for (x,y) = ({1},{2}) : {3}", new Object[]{this.getProjection().getClass().getName(), x, y, message});            
        } else {
            getProjection().getLogger().log(Level.FINE, "{0} - Solution not defined for (phi,theta) = ({1},{2}) : {3}", new Object[]{this.getProjection().getClass().getName(), x, y, message});            
        }
    }

    /**
     * Returns the X coordinate of the pixel.
     *
     * @return the x
     */
    public double getX() {
        return x;
    }

    /**
     * Returns the Y coordinate of the pixel.
     *
     * @return the y
     */
    public double getY() {
        return y;
    }

    /**
     * Returns true when the error is related to plane coordinate.
     * @return the isPlaneCoordinate
     */
    public boolean isIsPlaneCoordinate() {
        return isPlaneCoordinate;
    }

    @Override
    public String toString() {
        final String coordinates = isPlaneCoordinate ? "(x,y)" : "(phi,theta)";
        final String result;
        if (getMessage().isEmpty()) {
            result = MessageFormat.format("{0} - Solution not defined for {1} = ({2},{3})", this.getProjection().getClass().getName(), coordinates, getX(), getY());
        } else {
            result = MessageFormat.format("{0} - Solution not defined for {1} = ({2},{3}) : {4}", this.getProjection().getClass().getName(), coordinates, getX(), getY(), getMessage());
        }
        return result;
    }

}
