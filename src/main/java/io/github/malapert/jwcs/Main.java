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

import gnu.getopt.Getopt;
import gnu.getopt.LongOpt;
import io.github.malapert.jwcs.crs.Ecliptic;
import io.github.malapert.jwcs.crs.Equatorial;
import io.github.malapert.jwcs.datum.FK4;
import io.github.malapert.jwcs.datum.FK4NoEterms;
import io.github.malapert.jwcs.datum.FK5;
import io.github.malapert.jwcs.datum.ICRS;
import io.github.malapert.jwcs.datum.J2000;
import io.github.malapert.jwcs.position.SkyPosition;
import io.github.malapert.jwcs.crs.AbstractCrs;
import io.github.malapert.jwcs.crs.CrsFactory;
import io.github.malapert.jwcs.crs.gui.ConvertSelectionPanel;
import io.github.malapert.jwcs.proj.exception.JWcsException;
import io.github.malapert.jwcs.proj.exception.ProjectionException;
import io.github.malapert.jwcs.proj.gui.ProjectionSelectionPanel;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import nom.tam.fits.Fits;
import nom.tam.fits.FitsException;
import io.github.malapert.jwcs.datum.CoordinateReferenceFrame;

/**
 * The main class.
 *
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 */
public class Main {

    /**
     * Logger.
     */
    private final static Logger LOG = Logger.getLogger(Main.class.getName());

    /**
     * Default extension for FITS file.
     */
    private final static int DEFAULT_EXTENSION = 0;

    /**
     * List of programs.
     */
    private enum PROG {

        SKY_CONVERTER(""),
        PROJECT(""),
        UNPROJECT(""),
        GUI(null);

        private String commandLine;

        PROG(final String commandLine) {
            this.commandLine = commandLine;
        }

        public String getCommandLine() {
            return this.commandLine;
        }

        public void setCommandLine(final String commandLine) {
            this.commandLine = commandLine;
        }

    }

    /**
     * Exit code.
     */
    private enum EXIT {

        OK(0),
        USER_INPUT_ERROR(1),
        NO_EXIT(100);

        private final int code;

        EXIT(final int code) {
            this.code = code;
        }

        public int getCode() {
            return this.code;
        }
    }

    /**
     * Usage.
     */
    private static void usage() {
        final StringBuilder sb = new StringBuilder();
        sb.append("Usage: java -jar JWcs.jar -g PROG [OPTIONS]\n")
                .append("    or java -jar JWcs.jar --file HDR_FILE --project X,Y [OPTIONS]\n")
                .append("    or java -jar JWcs.jar --file HDR_FILE --unproject RA,DEC [OPTIONS]\n")
                .append("    or java -jar JWcs.jar --file HDR_FILE --convert RA,DEC --to SYS_TARGET [OPTIONS]\n")
                .append("    or java -jar JWcs.jar --convert RA,DEC --from SYS_ORGIN --to SYS_TARGET [OPTIONS]\n")
                .append("           where:\n")
                .append("               - PROG: either projection or converter\n")
                .append("               - HDR_FILE: Header FITS or FITS file\n")
                .append("               - X: pixel coordinate along X axis on the camera (starts to 1) \n")
                .append("               - Y: pixel coordinate along Y axis on the camera (starts to 1) \n")
                .append("               - RA: sky coordinate\n")
                .append("               - DEC: sky coordinate\n")
                .append("               - SYS_ORIGIN: sky system of the sky coordinates\n")
                .append("               - SYS_TARGET: convert sky coordinates to the SYS_TARGET\n\n")
                .append("           SYS_ORIGIN or SYS_TARGET can be:\n")
                .append("               - GALACTIC\n")
                .append("               - SUPER_GALACTIC\n")
                .append("               - EQUATORIAL\n")
                .append("               - EQUATORIAL(ICRS())\n")
                .append("               - EQUATORIAL(J2000())\n")
                .append("               - EQUATORIAL(FK5())\n")
                .append("               - EQUATORIAL(FK5(<equinox>))\n")
                .append("               - EQUATORIAL(FK4())\n")
                .append("               - EQUATORIAL(FK4(<equinox>))\n")
                .append("               - EQUATORIAL(FK4(<equinox>,<epoch>))\n")
                .append("               - EQUATORIAL(FK4_NO_E())\n")
                .append("               - EQUATORIAL(FK4_NO_E(<equinox>))\n")
                .append("               - EQUATORIAL(FK4_NO_E(<equinox>,<epoch>))\n")
                .append("               - ECLIPTIC\n")
                .append("               - ECLIPTIC(ICRS())\n")
                .append("               - ECLIPTIC(J2000())\n")
                .append("               - ECLIPTIC(FK5())\n")
                .append("               - ECLIPTIC(FK5(<equinox>))\n")
                .append("               - ECLIPTIC(FK4())\n")
                .append("               - ECLIPTIC(FK4(<equinox>))\n")
                .append("               - ECLIPTIC(FK4(<equinox>,<epoch>))\n")
                .append("               - ECLIPTIC(FK4_NO_E())\n")
                .append("               - ECLIPTIC(FK4_NO_E(<equinox>))\n")
                .append("               - ECLIPTIC(FK4_NO_E(<equinox>,<epoch>))\n")
                .append("\nProjection and sky conversion library\n")
                .append("\n")
                .append("Mandatory arguments to long options are mandatory for short options too.\n")
                .append("  -p, --project            Project a pixel to the sky\n")
                .append("  -u, --unproject          Unproject a point on the sky to 2D \n")
                .append("  -f, --file               Header file or Fits file starting by a scheme (ex: file://, http://)\n")
                .append("  -s, --from               Origin sky system\n")
                .append("  -t, --to                 Target sky system\n")
                .append("  -c, --convert            Convert a sky coordinate from a sky system to antoher one\n")
                .append("  -g, --gui                Display projection or converter with a GUI\n")
                .append("  -h, --help               Display this help and exit\n")
                .append("\n")
                .append("OPTIONS are the following:\n")
                .append("  -d, --debug              Sets the DEBUG level : ALL,CONFIG,FINER,FINEST,INFO,OFF,SEVERE,WARNING\n")
                .append("  -e, --extension          HDU number starting at 0 when --file argument is used. If not set, 0 is default\n")
                .append("  -r, --precision          Precision such as %.6f. By default, precision is set to %.15f\n");

        System.out.println(sb.toString());
        System.exit(EXIT.OK.getCode());
    }

    /**
     * Project camera to sky from command line.
     *
     * @param pos position in the camera frame
     * @param file file
     * @param extension HDU number when file is a FITS file
     * @param precision precision such as %.15f
     * @throws ProjectionException an error during the projection
     * @throws JWcsException JWcs error
     * @throws IOException an error when loading the FITS file
     * @throws URISyntaxException an error when loading the FITS file
     */
    private static void projectToSkyFromCommandLine(final String pos, final String file, final int extension, final String precision) throws ProjectionException, JWcsException, IOException, URISyntaxException {
        final String[] argumentsPos = pos.split(",");
        if (argumentsPos.length != 2)
            throw new IllegalArgumentException("The position " + pos + " is not correct");
        if (file == null)
            throw new IllegalArgumentException("--file argument is required");

        final URI uri = new URI(file);
        try {
            JWcsFits jf = new JWcsFits(new Fits(uri.toURL()), extension);
            jf.doInit();
            LOG.log(Level.INFO, "Executing pix2wcs(%s,%s)", argumentsPos);
            final double[] result = jf.pix2wcs(Double.valueOf(argumentsPos[0]), Double.valueOf(argumentsPos[1]));

            System.out.printf("(ra,dec)=(" + precision + ", " + precision + ")\n", result[0], result[1]);
            LOG.log(Level.INFO, "(ra,dec) = (%s,%s)", result);
        } catch (FitsException e) {
            throw new IOException(e.getMessage(), e);
        }
        
        /*
        final Map<String, String> keyMap = new HashMap<>();
        try (Fits fits = new Fits(uri.toURL())) {
            final Header hdr = fits.getHDU(extension).getHeader();
            final Cursor<String, HeaderCard> c = hdr.iterator();
            while (c.hasNext()) {
                final HeaderCard card = (HeaderCard) c.next();
                keyMap.put(card.getKey(), card.getValue());
            }
        } catch (nom.tam.fits.FitsException | IOException ex) {
            final HeaderFitsReader hdr = new HeaderFitsReader(uri.toURL());
            final List<List<String>> listKeywords = hdr.readKeywords();
            listKeywords.stream().forEach((keywordLine) -> {
                keyMap.put(keywordLine.get(0), keywordLine.get(1));
            });
        }

        final JWcsMap wcs = new JWcsMap(keyMap);
        wcs.doInit();
        LOG.log(Level.INFO, "Executing pix2wcs(%s,%s)", argumentsPos);
        final double[] result = wcs.pix2wcs(Double.valueOf(argumentsPos[0]), Double.valueOf(argumentsPos[1]));

        System.out.printf("(ra,dec)=(" + precision + ", " + precision + ")\n", result[0], result[1]);
        LOG.log(Level.INFO, "(ra,dec) = (%s,%s)", result);
        */
    }

    /**
     * Project sky to camera from command line.
     *
     * @param pos Sky position
     * @param file Header file
     * @param extension HDU number when file is a FITS file
     * @param precision precision such as %.15f
     * @throws JWcsException Exception
     * @throws IOException Exception
     * @throws URISyntaxException Exception
     */ 
    private static void projectToCameraFromCommandLine(final String pos, final String file, final int extension, final String precision) throws JWcsException, IOException, URISyntaxException {
        final String[] argumentsPos = pos.split(",");
        if (argumentsPos.length != 2)
            throw new IllegalArgumentException("The position " + pos + " is not correct");
        if (file == null)
            throw new IllegalArgumentException("--file argument is required");
        
        final URI uri = new URI(file);
        try {
            JWcsFits jf = new JWcsFits(new Fits(uri.toURL()), extension);
            jf.doInit();
            LOG.log(Level.INFO, "Executing wcs2pix(%s,%s)", argumentsPos);
            final double[] result = jf.wcs2pix(Double.valueOf(argumentsPos[0]), Double.valueOf(argumentsPos[1]));
            System.out.printf("(x,y)=(" + precision + ", " + precision + ")\n", result[0], result[1]);
            LOG.log(Level.INFO, "(x,y) = (%s,%s)", result);
        } catch (FitsException e) {
            throw new IOException(e.getMessage(), e);
        }

        /*
        final Map<String, String> keyMap = new HashMap<>();
        try (Fits fits = new Fits(uri.toURL())) {
            final Header hdr = fits.getHDU(extension).getHeader();
            final Cursor<String, HeaderCard> c = hdr.iterator();
            while (c.hasNext()) {
                final HeaderCard card = (HeaderCard) c.next();
                keyMap.put(card.getKey(), card.getValue());
            }
        } catch (nom.tam.fits.FitsException | IOException ex) {
            final HeaderFitsReader hdr = new HeaderFitsReader(uri.toURL());
            final List<List<String>> listKeywords = hdr.readKeywords();
            listKeywords.stream().forEach((keywordLine) -> {
                keyMap.put(keywordLine.get(0), keywordLine.get(1));
            });
        }
        final JWcsMap wcs = new JWcsMap(keyMap);
        wcs.doInit();
        LOG.log(Level.INFO, "Executing wcs2pix(%s,%s)", argumentsPos);
        final double[] result = wcs.wcs2pix(Double.valueOf(argumentsPos[0]), Double.valueOf(argumentsPos[1]));
        System.out.printf("(x,y)=(" + precision + ", " + precision + ")\n", result[0], result[1]);
        LOG.log(Level.INFO, "(x,y) = (%s,%s)", result);
        */
    }

    /**
     * Checks if the sky system has a reference system.
     *
     * @param skySystem sky system
     * @return True when the sky system has a reference frame
     */
    private static boolean hasReferenceFrame(final String skySystem) {
        return skySystem.contains("(");
    }

    /**
     * Extract reference frame name from coordinate reference system
     *
     * @param crs coordinate reference system
     * @return Reference frame name
     */
    private static String extractReferenceFrame(final String crs) {
        return crs.substring(crs.indexOf('(') + 1, crs.lastIndexOf(')'));
    }

    /**
     * Extracts parameters from the reference frame.
     *
     * @param referenceFrame reference frame
     * @return parameters of the reference frame
     */
    private static String[] extractReferenceFrameParameter(final String referenceFrame) {
        final String parameters = referenceFrame.substring(referenceFrame.indexOf('(') + 1, referenceFrame.lastIndexOf(')'));
        return parameters.isEmpty() ? null : parameters.split(",");
    }

    /**
     * Gets the reference frame with parameters
     *
     * @param refFrame reference frame type
     * @param parameters parameters of the reference frame
     * @return the reference frame
     */
    private static CoordinateReferenceFrame createReferenceFrameFactory(final CoordinateReferenceFrame.ReferenceFrame refFrame, final String[] parameters) {
        final CoordinateReferenceFrame ref;
        switch (refFrame) {
            case ICRS:
                ref = new ICRS();
                break;
            case FK5:
                ref = createFK5Frame(parameters);
                break;
            case FK4:
                ref = createFK4Frame(parameters);
                break;
            case FK4_NO_E:
                ref = createFK4NOEFrame(parameters);
                break;
            case J2000:
                ref = new J2000();
                break;
            default:
                throw new IllegalArgumentException("Reference frame "+refFrame+" not supported");
        }
        return ref;
    }

    /**
     * Creates a FK4 No eterms reference frame.
     * @param parameters parameters of the reference frame
     * @return FK4 No eterms reference frame
     */
    private static CoordinateReferenceFrame createFK4NOEFrame(final String[] parameters) {
        final CoordinateReferenceFrame ref;
        if (parameters != null) {
            switch (parameters.length) {
                case 1:
                    ref = new FK4NoEterms(parameters[0]);
                    break;
                case 2:
                    ref = new FK4NoEterms(parameters[0], parameters[1]);
                    break;
                default:
                    throw new IllegalArgumentException("");
            }
        } else {
            ref = new FK4NoEterms();
        }
        return ref;
    }

    /**
     * Creates a FK4 reference frame.
     * @param parameters parameters of the reference frame
     * @return FK4 reference frame
     */    
    private static CoordinateReferenceFrame createFK4Frame(final String[] parameters) {
        final CoordinateReferenceFrame ref;
        if (parameters != null) {
            switch (parameters.length) {
                case 1:
                    ref = new FK4(parameters[0]);
                    break;
                case 2:
                    ref = new FK4(parameters[0], parameters[1]);
                    break;
                default:
                    throw new IllegalArgumentException("");
            }
        } else {
            ref = new FK4();
        }
        return ref;
    }

    /**
     * Creates a FK5 reference frame.
     * @param parameters parameters of the reference frame
     * @return FK5 reference frame
     */    
    private static CoordinateReferenceFrame createFK5Frame(final String[] parameters) {
        final CoordinateReferenceFrame ref;
        if (parameters == null) {
            ref = new FK5();
        } else {
            ref = new FK5(parameters[0]);
        }
        return ref;
    }

    /**
     * Gets the coordinate reference system.
     *
     * @param crs the coordinate reference system name
     * @return the AbstractCrs object
     */
    private static AbstractCrs getCrs(final String crs) {
        final AbstractCrs result;
        if (hasReferenceFrame(crs)) {
            final String skySystemName = crs.substring(0, crs.indexOf('('));
            final String refFrameStr = extractReferenceFrame(crs);
            final String refFrameName = refFrameStr.substring(0, refFrameStr.indexOf('('));
            final CoordinateReferenceFrame.ReferenceFrame refFrame = CoordinateReferenceFrame.ReferenceFrame.valueOf(refFrameName);
            final String[] parameters = extractReferenceFrameParameter(refFrameStr);
            final CoordinateReferenceFrame refSystem = createReferenceFrameFactory(refFrame, parameters);
            result = CrsFactory.create(AbstractCrs.CoordinateReferenceSystem.valueOf(skySystemName));
            switch (result.getCoordinateReferenceSystem()) {
                case EQUATORIAL:
                    ((Equatorial) result).setCoordinateReferenceFrame(refSystem);
                    break;
                case ECLIPTIC:
                    ((Ecliptic) result).setCoordinateReferenceFrame(refSystem);
                    break;
                default:
                    throw new IllegalArgumentException("Unknown sky system");
            }
        } else {
            result = CrsFactory.create(AbstractCrs.CoordinateReferenceSystem.valueOf(crs));
        }
        return result;
    }

    /**
     * Converts a position frome a coordinate reference system to another one.
     *
     * @param pos Sky position
     * @param file Fits or header file
     * @param from source crs
     * @param to target crs
     * @param extension FITS extension
     * @param precision precision such as %.15f
     * @throws URISyntaxException Cannot retrieve the Header file
     * @throws IOException Header file not found
     * @throws JWcsException JWS Error
     * @throws IllegalArgumentException Either --file argument or --from and
     * --to arguments are required
     */
    private static void convertFromCommandLine(final String pos, final String file, final String from, final String to, final int extension, final String precision) throws URISyntaxException, IOException, JWcsException {
        final AbstractCrs crsFrom;
        if (file == null && from == null && to == null)
            throw new IllegalArgumentException("Either --file argument or --from and --to arguments are required");
        
        if (file != null) {
            final URI uri = new URI(file);
            try {
                JWcsFits jf = new JWcsFits(new Fits(uri.toURL()), extension);
                jf.doInit();
                crsFrom = jf.getCrs();
            } catch (FitsException e) {
                throw new IOException(e.getMessage(), e);
            }            
            
            /*
            final Map<String, String> keyMap = new HashMap<>();
            try {
                final Fits fits = new Fits(uri.toURL());
                final Header hdr = fits.getHDU(extension).getHeader();
                final Cursor<String, HeaderCard> c = hdr.iterator();
                while (c.hasNext()) {
                    final HeaderCard card = (HeaderCard) c.next();
                    keyMap.put(card.getKey(), card.getValue());
                }
            } catch (nom.tam.fits.FitsException | IOException ex) {
                final HeaderFitsReader hdr = new HeaderFitsReader(uri.toURL());
                final List<List<String>> listKeywords = hdr.readKeywords();
                listKeywords.stream().forEach((keywordLine) -> {
                    keyMap.put(keywordLine.get(0), keywordLine.get(1));
                });
            }
            final JWcsMap wcs = new JWcsMap(keyMap);
            wcs.doInit();
            crsFrom = wcs.getCrs();
            */
        } else {
            crsFrom = getCrs(from);
        }
        
        final String crsTarget = to;
        final double[] skyPos = Arrays.stream(pos.split(","))
                .mapToDouble(Double::parseDouble)
                .toArray();
        final AbstractCrs crsTo = getCrs(crsTarget);
        LOG.log(Level.INFO, "Executing convertTo(%s, %s,%s)", new Object[]{crsTo, skyPos[0], skyPos[1]});
        final SkyPosition skyPosition = crsFrom.convertTo(crsTo, skyPos[0], skyPos[1]);
        System.out.printf(precision + ", " + precision + "\n", skyPosition.getLongitude(), skyPosition.getLatitude());
        LOG.log(Level.INFO, "(longitude,latitude) = (%s,%s)", skyPosition);
    }

    /**
     * Sets debug level.
     *
     * @param rootLogger the root logger
     * @param levelArg the debug level
     */
    private static void setDebug(final Logger rootLogger, final Getopt levelArg) {
        final Level levelInfo = Level.parse(levelArg.getOptarg());
        final Handler[] handlers = rootLogger.getHandlers();
        rootLogger.setLevel(levelInfo);
        for (final Handler h : handlers) {
            if (h instanceof FileHandler) {
                h.setLevel(levelInfo);
            } else if (h instanceof ConsoleHandler) {
                h.setLevel(levelInfo);
            }
        }
        rootLogger.setLevel(levelInfo);
    }

    /**
     * The program starts with this method.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        EXIT returnedCode = EXIT.OK;
        boolean isGui = false;
        int c;
        int extension = DEFAULT_EXTENSION;
        String from = null;
        String to = null;
        String file = null;
        String precision = "%.15f";
        String progGui = null;
        final List<PROG> progChoice = new ArrayList<>();
        LongOpt[] longopts = new LongOpt[11];
        final Logger rootLogger = Logger.getLogger("");
        rootLogger.setLevel(Level.OFF);

        longopts[0] = new LongOpt("help", LongOpt.NO_ARGUMENT, null, 'h');
        longopts[1] = new LongOpt("gui", LongOpt.REQUIRED_ARGUMENT, null, 'g');
        longopts[2] = new LongOpt("project", LongOpt.REQUIRED_ARGUMENT, null, 'p');
        longopts[3] = new LongOpt("unproject", LongOpt.REQUIRED_ARGUMENT, null, 'u');
        longopts[4] = new LongOpt("convert", LongOpt.REQUIRED_ARGUMENT, null, 'c');
        longopts[5] = new LongOpt("debug", LongOpt.REQUIRED_ARGUMENT, null, 'd');
        longopts[6] = new LongOpt("file", LongOpt.REQUIRED_ARGUMENT, null, 'f');
        longopts[7] = new LongOpt("from", LongOpt.REQUIRED_ARGUMENT, null, 's');
        longopts[8] = new LongOpt("to", LongOpt.REQUIRED_ARGUMENT, null, 't');
        longopts[9] = new LongOpt("extension", LongOpt.REQUIRED_ARGUMENT, null, 'e');
        longopts[10] = new LongOpt("precision", LongOpt.REQUIRED_ARGUMENT, null, 'r');
        // 
        final Getopt g = new Getopt("JWcs", args, "-::p:u:c:d:f:s:t:e:r:g:h;", longopts);
        g.setOpterr(true);
        //
        while ((c = g.getopt()) != -1) {
            switch (c) {
                case 'r':
                    precision = g.getOptarg();
                    break;
                case 'e':
                    try {
                        extension = Integer.parseInt(g.getOptarg());
                    } catch (NumberFormatException ex) {
                        System.err.println(ex);
                        returnedCode = EXIT.USER_INPUT_ERROR;
                    }
                    break;
                case 'f':
                    file = g.getOptarg();
                    break;
                case 's':
                    from = g.getOptarg();
                    break;
                case 't':
                    to = g.getOptarg();
                    break;
                case 'g':
                    progChoice.add(PROG.GUI);
                    progGui = g.getOptarg();
                    break;
                case 'c':
                    PROG.SKY_CONVERTER.setCommandLine(g.getOptarg());
                    progChoice.add(PROG.SKY_CONVERTER);
                    break;
                case 'd':
                    setDebug(rootLogger, g);
                    break;
                case 'p':
                    PROG.PROJECT.setCommandLine(g.getOptarg());
                    progChoice.add(PROG.PROJECT);
                    break;
                case 'u':
                    PROG.UNPROJECT.setCommandLine(g.getOptarg());
                    progChoice.add(PROG.UNPROJECT);
                    break;
                case 'h':
                    usage();
                    break;
                case ':':
                    System.err.println("You need an argument for option "
                            + (char) g.getOptopt());
                    LOG.log(Level.SEVERE, "You need an argument for option %s", g.getOptopt());
                    returnedCode = EXIT.USER_INPUT_ERROR;
                    break;
                //
                case '?':
                    System.err.println("The option '" + (char) g.getOptopt()
                            + "' is not valid");
                    LOG.log(Level.SEVERE, "The option %s is not valid", g.getOptopt());
                    returnedCode = EXIT.USER_INPUT_ERROR;
                    break;
                //
                default:
                    System.err.println("getopt() returned " + c);
                    LOG.log(Level.SEVERE, "getopt() returned %s", c);
                    returnedCode = EXIT.USER_INPUT_ERROR;
                    break;
            }
        }
        for (int i = g.getOptind(); i < args.length; i++) {
            System.err.println("Non option argv element: " + args[i] + "\n");
            LOG.log(Level.SEVERE, "Non option argv element: %s", args[i]);
            returnedCode = EXIT.USER_INPUT_ERROR;
        }
        if (EXIT.OK != returnedCode) {
            System.exit(returnedCode.getCode());
        }

        if (progChoice.size() != 1) {
            System.err.println("You need to select only one of the available program : gui, project, unproject, converter");
            returnedCode = EXIT.USER_INPUT_ERROR;
            usage();
            System.exit(returnedCode.getCode());
        }

        try {

            final PROG prog = progChoice.get(0);
            switch (prog) {
                case GUI:
                    if (null != progGui) {
                        switch (progGui) {
                            case "projection":
                                java.awt.EventQueue.invokeLater(() -> {
                                    try {
                                        ProjectionSelectionPanel.createWindow();
                                    } catch (IOException ex) {
                                        System.err.println("Error : " + ex.getMessage());
                                        LOG.log(Level.SEVERE, null, ex);
                                    }
                                });
                                break;
                            case "converter":
                                java.awt.EventQueue.invokeLater(() -> {
                                    ConvertSelectionPanel.createWindow();
                                });
                                break;
                            default:
                                throw new IllegalArgumentException("The GUI program " + progGui + " is not supported");
                        }
                    }
                    isGui = true;
                    break;
                case PROJECT:
                    projectToSkyFromCommandLine(prog.getCommandLine(), file, extension, precision);
                    break;
                case UNPROJECT:
                    projectToCameraFromCommandLine(prog.getCommandLine(), file, extension, precision);
                    break;
                case SKY_CONVERTER:
                    convertFromCommandLine(prog.getCommandLine(), file, from, to, extension, precision);
                    break;
                default:
                    throw new IllegalArgumentException(prog.name() + " not supported");
            }
            returnedCode = EXIT.OK;
        } catch (JWcsException | IOException | RuntimeException | URISyntaxException ex) {
            LOG.log(Level.SEVERE, "Error : ", ex.getMessage());
            System.err.println("Error: " + ex.getMessage());
            returnedCode = EXIT.USER_INPUT_ERROR;
        } finally {
            if (!isGui) {
                System.exit(returnedCode.getCode());
            }
        }
    }
}
