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
package io.github.malapert.jwcs;

import io.github.malapert.jwcs.proj.exception.JWcsError;
import io.github.malapert.jwcs.proj.exception.JWcsException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.util.FastMath;

/**
 * A class allowing to compute the WCS by initializing the WCS structure by
 * reading a map containing the WCS keywords.
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class JWcsMap extends AbstractJWcs {

    /**
     * WCS keywords provided by the user.
     */
    private Map keywords;

    /**
     * Clones a JWcsMap.
     *
     * @param wcs JWcsMap object
     */
    public JWcsMap(final JWcsMap wcs) {
        setKeywords(wcs.keywords);
    }

    /**
     * Creates a AbstractJWcs bases on a map.
     *
     * @param keywords map of keywords
     */
    public JWcsMap(final Map keywords) {
        this.keywords = keywords;
    }

    @Override
    public void doInit() throws JWcsException {
        super.init();
    }

    @Override
    public int wcsaxes() {
        return this.getValueAsInt("NAXIS");
    }

    @Override
    public boolean hasKeyword(final String keyword) {
        return this.getKeywords().containsKey(keyword);
    }

    /**
     * Get the <CODE>int</CODE> value associated with the given key.
     *
     * @param keyword The header key.
     * @return The associated value.
     * @throws JWcsError when the keyword is not found
     */
    @Override
    public int getValueAsInt(final String keyword) {
        final String result = getValueAsString(keyword);
        if (result == null) {
            throw new JWcsError(keyword + " not found");
        } else {
            return Integer.parseInt(result);
        }
    }

    /**
     * Get the <CODE>double</CODE> value associated with the given key.
     *
     * @param keyword The header key.
     * @return The associated value.
     * @throws JWcsError when the keyword is not found
     */
    @Override
    public double getValueAsDouble(final String keyword) {
        final String result = getValueAsString(keyword);
        if (result == null) {
            throw new JWcsError(keyword + " not found");
        } else {
            return Double.valueOf(result);
        }
    }

    /**
     * Get the <CODE>float</CODE> value associated with the given key.
     *
     * @param keyword The header key.
     * @return The associated value.
     * @throws JWcsError when the keyword is not found
     */
    @Override
    public float getValueAsFloat(final String keyword) {
        final String result = getValueAsString(keyword);
        if (result == null) {
            throw new JWcsError(keyword + " not found");
        } else {
            return Float.valueOf(result);
        }
    }

    /**
     * Get the <CODE>String</CODE> value associated with the given key.
     *
     * @param keyword The header key.
     * @return The associated value or null when the keyword is not found
     */
    @Override
    public String getValueAsString(final String keyword) {
        return (String) this.getKeywords().getOrDefault(keyword, null);
    }

    @Override
    public Iterator iterator() {
        return this.getKeywords().keySet().iterator();
    }

    /**
     * Returns all the keywords.
     *
     * @return the keywords
     */
    public Map getKeywords() {
        return keywords;
    }

    /**
     * Sets the whole map.
     *
     * @param keywords the keywords to set
     */
    private void setKeywords(final Map keywords) {
        this.keywords = keywords;
    }

    /**
     * Gets a projection based on a projection code.
     *
     * @param projectionCode projection code
     * @return Returns a AbstractJWcs projection
     * @throws JWcsException When a JWcsException happens
     */
    public static AbstractJWcs getProjection(final String projectionCode) throws JWcsException {
        final Map wcsKeywords = new HashMap();
        wcsKeywords.put(AbstractJWcs.NAXIS, "2");
        wcsKeywords.put(AbstractJWcs.NAXIS1, "600");
        wcsKeywords.put(AbstractJWcs.NAXIS2, "300");
        wcsKeywords.put(AbstractJWcs.RADESYS, "ICRS");
        wcsKeywords.put(AbstractJWcs.CRPIX1, "300");
        wcsKeywords.put(AbstractJWcs.CRPIX2, "150");
        wcsKeywords.put(AbstractJWcs.CRVAL1, "0");
        if (null != projectionCode) {
            fillProjectionParameters(projectionCode, wcsKeywords);
        }
        wcsKeywords.put(AbstractJWcs.CD11, String.valueOf(180d / 300d));
        wcsKeywords.put(AbstractJWcs.CD12, "0");
        wcsKeywords.put(AbstractJWcs.CD21, "0");
        wcsKeywords.put(AbstractJWcs.CD22, String.valueOf(90d / 150d));
        wcsKeywords.put(AbstractJWcs.CTYPE1, "RA---" + projectionCode);
        wcsKeywords.put(AbstractJWcs.CTYPE2, "DEC--" + projectionCode);
        final AbstractJWcs wcs = new JWcsMap(wcsKeywords);
        wcs.doInit();
        return wcs;
    }
    
    /**
     * Fills projection parameters.
     * 
     * @param projectionCode projection code
     * @param wcsKeywords WCS keywords
     */
    private static void fillProjectionParameters(final String projectionCode, final Map wcsKeywords) {
        switch (projectionCode) {
            case "ZEA":
            case "TAN":
            case "TPV":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "90");
                break;
            case "BON":
                wcsKeywords.put(AbstractJWcs.PV21, "45");
                wcsKeywords.put(AbstractJWcs.CRVAL2, "0");
                break;
            case "SZP":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "0");
                wcsKeywords.put(AbstractJWcs.PV21, "2");
                wcsKeywords.put(AbstractJWcs.PV22, "180");
                wcsKeywords.put(AbstractJWcs.PV23, "60");
                break;
            case "STG":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "-90");
                break;
            case "CYP":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "0");
                wcsKeywords.put(AbstractJWcs.PV21, "1");
                wcsKeywords.put(AbstractJWcs.PV22, String.valueOf(FastMath.sqrt(2) * 0.5));
                break;
            case "COP":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "90");
                wcsKeywords.put(AbstractJWcs.PV21, "45");
                wcsKeywords.put(AbstractJWcs.PV22, "25");
                break;
            case "COO":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "90");
                wcsKeywords.put(AbstractJWcs.PV21, "45");
                wcsKeywords.put(AbstractJWcs.PV22, "25");
                break;
            case "COE":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "-90");
                wcsKeywords.put(AbstractJWcs.PV21, "-45");
                wcsKeywords.put(AbstractJWcs.PV22, "25");
                break;
            case "COD":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "90");
                wcsKeywords.put(AbstractJWcs.PV21, "45");
                wcsKeywords.put(AbstractJWcs.PV22, "25");
                break;
            case "AZP":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "60");
                wcsKeywords.put(AbstractJWcs.PV21, "2");
                wcsKeywords.put(AbstractJWcs.PV22, "30");
                break;
            case "ARC":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "90");
                break;
            case "ZPN":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "90");
                wcsKeywords.put(AbstractJWcs.PV20, "0.050");
                wcsKeywords.put(AbstractJWcs.PV21, "0.975");
                wcsKeywords.put(AbstractJWcs.PV22, "-0.807");
                wcsKeywords.put("PV2_3", "0.337");
                wcsKeywords.put("PV2_4", "-0.065");
                wcsKeywords.put("PV2_5", "0.010");
                wcsKeywords.put("PV2_6", "0.003");
                wcsKeywords.put("PV2_7", "-0.001");
                break;
            case "AIR":
                wcsKeywords.put(AbstractJWcs.CRVAL2, "90");
                wcsKeywords.put(AbstractJWcs.PV21, "45");
                break;
            default:
                wcsKeywords.put(AbstractJWcs.CRVAL2, "0");
                break;
        }     
    }

    /**
     * Returns True when all required WCS keywords are there otherwise False.
     *
     * @return True when all required WCS keywords are there otherwise False
     */
    private boolean hasRequiredWcs() {
        final boolean hasOldCd = this.hasKeyword(CDELT1) && this.hasKeyword(CDELT2) && this.hasKeyword(CROTA2);
        final boolean hasCd = hasCd() || hasOldCd;
        return this.hasKeyword(CTYPE1) && this.hasKeyword(CTYPE2)
                && this.hasKeyword(CRPIX1) && this.hasKeyword(CRPIX2)
                && this.hasKeyword(CRVAL1) && this.hasKeyword(CRVAL2)
                && hasCd;
    }

    @Override
    protected void checkWcs() throws JWcsException {
        if(!hasRequiredWcs()) {
            fillStackError();
        }
    }
    
    /**
     * Fills stack of error in JWcsException.
     * @throws JWcsException All missing mandatory keywords
     */
    private void fillStackError() throws JWcsException {
        final List<String> check = new ArrayList();
        final boolean hasOldCd = this.hasKeyword(CDELT1) && this.hasKeyword(CDELT2) && this.hasKeyword(CROTA2);
        final boolean hasCd = hasCd() || hasOldCd;            
        if (!hasCd) {
            check.add("(CDELT1,CDELT2,CROTA2) or (CD11,C12,CD21,CD22) or (PC11, PC12, P21, PC22, CDELT1, CDELT2) are missing");
        }
        if (!this.hasKeyword(CTYPE1)) {
            check.add(CTYPE1 + " is missing");
        }
        if (!this.hasKeyword(CTYPE2)) {
            check.add(CTYPE2 + " is missing");
        }
        if (!this.hasKeyword(CRPIX1)) {
            check.add(CRPIX1 + " is missing");
        }
        if (!this.hasKeyword(CRPIX2)) {
            check.add(CRPIX2 + " is missing");
        }
        if (!this.hasKeyword(CRVAL1)) {
            check.add(CRVAL1 + " is missing");
        }
        if (!this.hasKeyword(CRVAL2)) {
            check.add(CRVAL2 + " is missing");
        }
        throw new JWcsException(check.toString());
        
    }
}
