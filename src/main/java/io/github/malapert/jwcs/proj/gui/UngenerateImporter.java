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

import java.io.*;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

/**
 * An importer for the ESRI Ungenerate file format.
 *
 * @author Bernhard Jenny, Institute of Cartography, ETH Zurich.
 */
public final class UngenerateImporter {

    private UngenerateImporter() {
    }

    /**
     * Reads an Ungenerate file and returns the found lines in a Vector.
     *
     * @param filePath The file to import.
     * @return A Vector containing all read MapLines.
     * @throws java.io.IOException Exception
     */
    public static List<MapLine> importData(final String filePath) throws IOException {
        return importData(new FileInputStream(filePath));
    }

    /**
     * Returns a list of lines from the import.
     *
     * @param inputStream input stream
     * @return a list of lines
     * @throws IOException File not found
     */
    public static List<MapLine> importData(final InputStream inputStream) throws IOException {
        // store the read lines in this vector.
        final List<MapLine> lines = new ArrayList<>();

        // read lines until we reach the end of the file
        try ( // open the file
              BufferedReader in = new BufferedReader(new InputStreamReader(inputStream, Charset.defaultCharset()))) {
            while (in.readLine() != null) {
                final MapLine line = readLine(in);
                lines.add(line);

            }
        }

        // return the read lines
        return lines;
    }

    /**
     * Reads a MapLine from the current file position.
     *
     * @param in The reader that provides the data to read.
     * @return MapLine The read line.
     * @exception java.io.IOException Exception
     */
    private static MapLine readLine(final BufferedReader in) throws java.io.IOException {
        String str;
        final MapLine line = new MapLine();

        while (true) {
            str = in.readLine();
            if (str == null || str.isEmpty()) {
                break;
            }
            str = str.trim().toLowerCase();
            if (str.startsWith("end")) {
                break;
            }
            try {
                final StringTokenizer tokenizer = new StringTokenizer(str, " \t,");
                double x = Double.parseDouble(tokenizer.nextToken());
                final double y = Double.parseDouble(tokenizer.nextToken());
                if (x < 0) {
                    x = 360 + x;
                }
                line.addPoint(x, y);
            } catch (NoSuchElementException e) {
                // found a line without any readable data. Just read the next line
            }
        }
        return line;
    }

}
