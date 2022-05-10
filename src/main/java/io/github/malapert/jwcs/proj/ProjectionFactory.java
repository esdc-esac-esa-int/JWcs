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

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 */
public class ProjectionFactory {
    
    /**
     * LOG.
     */
    protected final static Logger LOG = Logger.getLogger(ProjectionFactory.class.getName());
    
    
    private static Map registry = new HashMap();
    private static Map nameMap = new HashMap();    
    
    public static AbstractProjection getNamedProjection(String name) {    
        if (registry.isEmpty()) {
            initialize();
        }
        String projectionCode= (String) nameMap.get(name);
        return getNamedCodeProjection(projectionCode);        
    }
    
    public static Object[] getOrderedProjectionCodes() {
        if (registry.isEmpty()) {
            initialize();
        }
        Object[] names = nameMap.keySet().toArray();
        Arrays.sort(names);
        return names;        
    }
    
    public static AbstractProjection getNamedCodeProjection(String name) {
        if (registry.isEmpty()) {
            initialize();
        }
        Class cls = (Class) registry.get(name);
        if (cls != null) {
            try {
                AbstractProjection projection = (AbstractProjection) cls.newInstance();
                return projection;
            } catch (IllegalAccessException | InstantiationException e) {
            }
        }
        return null;        
    }
    
    private static void initialize() {
        try {
            register("AIR", AIR.class);
            register("AIT", AIT.class);
            register("ARC", ARC.class);
            register("AZP", AZP.class);
            register("BON", BON.class);
            register("CAR", CAR.class);
            register("CEA", CEA.class);
            register("COD", COD.class);
            register("COE", COE.class);
            register("COO", COO.class);
            register("COP", COP.class);
            register("CYP", CYP.class);            
            register("MER", MER.class);
            register("MOL", MOL.class);
            register("NCP", NCP.class);
            register("PAR", PAR.class);
            register("PCO", PCO.class);
            register("SFL", SFL.class);   
            register("SIN", SIN.class);
            register("STG", STG.class);
            register("SZP", SZP.class);
            register("TAN", TAN.class);
            register("TPV", TPV.class);
            register("ZEA", ZEA.class);
            register("ZPN", ZPN.class);              
        } catch (NoSuchMethodException | InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException ex) {
            Logger.getLogger(ProjectionFactory.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    private static void register(String projectionCode, Class cls) throws NoSuchMethodException, InstantiationException, IllegalAccessException, IllegalArgumentException, InvocationTargetException {
        registry.put(projectionCode, cls);
        final AbstractProjection projection;
        final Constructor<?> constructor = cls.getConstructor();
        final Object instance = constructor.newInstance();
        projection = (AbstractProjection) instance;
        nameMap.put(projectionCode, projection.getName());   
    }
    
}
