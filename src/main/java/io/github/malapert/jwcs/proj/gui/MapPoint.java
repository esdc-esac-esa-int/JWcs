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
package io.github.malapert.jwcs.proj.gui;

/**
 *
 * @author malapert
 */
public final class MapPoint {

    /**
     * X coordinate along X axis.
     */
    public double x;
    
    /**
     * Y coordinate along Y axis.
     */
    public double y;
    
    /**
     * Constructor.
     */
    public MapPoint() {
        //do nothing
    }

    /**
     * Constructor.
     * @param x coordinate along X axis
     * @param y coordinate along Y axis
     */
    public MapPoint (final double x, final double y) {
        this.x = x;
        this.y = y;
    }

    
    @Override
    public boolean equals(final Object obj) {
        final boolean isEqual;
        if (obj != null && obj instanceof MapPoint p) {
            isEqual = x == p.x && y == p.y;
        } else {
            isEqual = false;
        } 
        return isEqual;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 89 * hash + (int) (Double.doubleToLongBits(this.x) ^ (Double.doubleToLongBits(this.x) >>> 32));
        hash = 89 * hash + (int) (Double.doubleToLongBits(this.y) ^ (Double.doubleToLongBits(this.y) >>> 32));
        return hash;
    }
}
