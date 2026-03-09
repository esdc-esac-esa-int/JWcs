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

package io.github.malapert.jwcs.utility;

import java.io.Serializable;
import java.text.NumberFormat;
import java.util.Locale;
import java.util.StringTokenizer;


/**
 * Class representing a value of the form "hours:min:sec".
 *
 * @author Allan Brighton
 * @version $Revision: 1.4 $
 */
public final class HMS implements Serializable {

    /**
     * On the handling of -0: from the javadoc for Double.equals():
     * "If d1 represents +0.0 while d2 represents -0.0, or vice versa,
     * the equal test has the value false, even though +0.0==-0.0 has the
     * value true."
     * The test for 0.0 != -0.0 only works with Double.equals(minusZero).
     *
     * <p>This case shows up in HMS values with zero hours and negative values,
     * such as "-00 24 32"
     */
    private static final Double MINUS_ZERO = -0.0;

    /**
     * Number formats for 2 digit hours and minutes.
     */
    private static final NumberFormat NF = NumberFormat.getInstance(Locale.US);

    /**
     * Number formats for seconds.
     */
    private static final NumberFormat NF_SEC = NumberFormat.getInstance(Locale.US);
    private static final long serialVersionUID = 6425466963081211760L;

    private static final double EPSILON = 1E-14;

    static {
        NF.setMinimumIntegerDigits(2);
        NF.setMaximumIntegerDigits(2);
        NF.setMaximumFractionDigits(0);

        NF_SEC.setMinimumIntegerDigits(2);
        NF_SEC.setMaximumIntegerDigits(2);
        NF_SEC.setMinimumFractionDigits(3);
        NF_SEC.setMaximumFractionDigits(3);
    }

    /**
     * number of hours.
     */
    private int hours;

    /**
     * number of minutes.
     */
    private int min;

    /**
     * number of seconds.
     */
    private double sec;

    /**
     * value converted to decimal.
     */
    private double val;

    /**
     * set to 1 or -1.
     */
    private byte sign = 1;

    /**
     * true if value has been initialized .
     */
    private boolean initialized;


    /**
     * Default constructor: initialize to null values.
     */
    public HMS() {
        this.initialized = false;
    }

    /**
     * Initialize with the given hours, minutes and seconds.
     * @param hours hours
     * @param min minutes
     * @param sec seconds
     */
    public HMS(final double hours, final int min, final double sec) {
        set(hours, min, sec);
    }

    /**
     * Initialize from a decimal hours value and calculate H:M:S.sss.
     * @param val the HMS as double
     */
    public HMS(final double val) {
        setVal(val);
    }

    /**
     * Copy constructor.
     * @param hms HDS as hms
     */
    public HMS(final HMS hms) {
        setVal(hms.val);
    }

    /**
     * Initialize from a string value, in format H:M:S.sss, hh, or H M
     * S. 
     * 
     * <p>If the value is not in H:M:S and is not an integer (has a
     * decimal point), assume the value is in deg convert to hours by
     * dividing by 15. (Reason: some catalog servers returns RA in h:m:value
     * while others return it in decimal deg.)
     * @param value HMS as String
     */
    public HMS(final String value) {
        this(value, false);
    }

    /**
     * Initialize from a string value, in format H:M:S.sss, hh, or
     * H M S.  
     * 
     * <p>If the value is not in H:M:S and is not an
     * integer (has a decimal point), and hflag is true,
     * assume the value is in deg and convert to hours by dividing by 15.
     *
     * @param value the RA string
     * @param hflag if true, assume RA is always in hours, otherwise, if it has a decimal point,
     * assume deg
     * @throws IllegalArgumentException Expected a string of the form hh:mm:ss.sss
     */
    public HMS(final String value, final boolean hflag) {
        final String valueProcessed = value.replace(",", "."); // Treat ',' like '.', by request
        double[] vals = {0.0, 0.0, 0.0};
        final StringTokenizer tok = new StringTokenizer(valueProcessed, ": ");
        int n = 0;
        while (n < 3 && tok.hasMoreTokens()) {
            vals[n++] = Double.parseDouble(tok.nextToken());
        }

        if (n >= 2) {
            set(vals[0], (int) vals[1], vals[2]);
        } else if (n == 1) {
            if (!hflag && value.indexOf('.') != -1) {
                setVal(vals[0] / 15.);
            } else {
                setVal(vals[0]);
            }
        } else {
            throw new IllegalArgumentException("Expected a string of the form hh:mm:ss.sss, but got: '" + value + "'");
        }
    }

    /**
     * Set the hours, minutes and seconds.
     * @param hours hours
     * @param min minutes
     * @param sec seconds
     */
    public void set(final double hours, final int min, final double sec) {
        this.hours = (int) hours;
        this.min = min;
        this.sec = sec;

        val = (sec / 60.0 + min) / 60.0;

        if (hours < 0.0 || Double.compare(hours, MINUS_ZERO) == 0) {
            val = hours - val;
            this.hours = -this.hours;
            sign = -1;
        } else {
            val = this.hours + val;
            sign = 1;
        }
        initialized = true;
    }

    /**
     * Set from a decimal value (hours) and calculate H:M:S.sss.
     * @param val value
     */
    public void setVal(final double val) {
        this.val = val;

        double v = val; // check also for neg zero
        if (v < 0.0 || Double.compare(v, MINUS_ZERO) == 0) {
            sign = -1;
            v = -v;
        } else {
            sign = 1;
        }

        final double dd = v + 0.0000000001;
        hours = (int) dd;
        final double md = (dd - hours) * 60.;
        min = (int) md;
        sec = (md - min) * 60.;
        initialized = true;
    }

    /**
     * Return the value as a String in the form hh:mm:ss.sss.
     * 
     * <p>Seconds are formatted with leading zero if needed.
     * The seconds are formatted with 3 digits of precision.
     */
    @Override
    public String toString() {
        final String secs = NF_SEC.format(sec);

        // sign
        final String signStr;
        if (sign == -1) {
            signStr = "-";
        } else {
            signStr = "";
        }

        return signStr
                + NF.format(hours)
                + ":"
                + NF.format(min)
                + ":"
                + secs;
    }

    /**
     * Return the value as a String in the form hh:mm:ss.sss,
     * or if showSeconds is false, hh:mm.
     * @param showSeconds show seconds
     * @return the value as a String in the form hh:mm[:ss.sss]
     */
    public String toString(final boolean showSeconds) {
        if (showSeconds) {
            return toString();
        }

        // sign
        final String signStr;
        if (sign == -1) {
            signStr = "-";
        } else {
            signStr = " ";
        }

        return signStr
                + NF.format(hours)
                + ":"
                + NF.format(min);
    }

    /**
     * Return true if this object has been initialized with a valid value.
     * @return true if this object has been initialized with a valid value
     */
    public boolean isInitialized() {
        return initialized;
    }

    /**
     * Return the number of hours (not including minutes or seconds).
     * @return the number of hours
     */
    public int getHours() {
        return hours;
    }

    /**
     * Return the number of minutes (not including hours or seconds).
     * @return the number of minutes
     */
    public int getMin() {
        return min;
    }

    /**
     * Return the number of seconds (not including hours and minutes).
     * @return he number of seconds
     */
    public double getSec() {
        return sec;
    }

    /**
     * Return the value (fractional number of hours) as a double.
     * @return the value as double
     */
    public double getVal() {
        return val;
    }

    /**
     * Return the sign of the value.
     * @return the sign of the value
     */
    public byte getSign() {
        return sign;
    }

    /**
     * Define equality based on the value.
     */
    @Override
    public int hashCode() {
        int hash = 7;
        hash = 97 * hash + (int) (Double.doubleToLongBits(this.val) ^ (Double.doubleToLongBits(this.val) >>> 32));
        return hash;
    }

    @Override
    public boolean equals(final Object obj) {
        return obj instanceof HMS && Math.abs(val - ((HMS) obj).val) <= EPSILON;
    }
}
