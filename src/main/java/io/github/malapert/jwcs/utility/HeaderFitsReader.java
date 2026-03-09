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

import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;

/**
 * Reads data from a HeaderFits map.
 * @author Jean-Christophe Malapert (jcmalapert@gmail.com)
 * @version 2.0
 */
public class HeaderFitsReader {

    /**
     * Separator between keyword/value.
     */
    private static final String SEPARATOR = "=";

    /**
     * Reader.
     */
    private final Reader source;

    /**
     * Constructs a HeaderFitsReader based on a reader.
     * @param source reader
     */
    public HeaderFitsReader(final Reader source) {
        this.source = source;
    }

    /**
     * Constructs a HeaderFitsReader based on a file.
     * @param file file
     * @throws FileNotFoundException File not found
     */
    public HeaderFitsReader(final File file) throws FileNotFoundException {
        this.source = new FileReader(file);
    }

    /**
     * Constructs a HeaderFitsReader based on a URI.
     * @param uri uri of the file
     * @throws FileNotFoundException  File not found
     */
    public HeaderFitsReader(final URI uri) throws FileNotFoundException {
        this.source = new FileReader(new File(uri));
    }

    /**
     * Constructs a HeaderFitsReader based on a URL.     
     * @param url URL
     * @throws IOException File not found
     */
    public HeaderFitsReader(final URL url) throws IOException {
        this.source = new InputStreamReader(url.openStream(), Charset.defaultCharset());
    }

    /**
     * Constructs a HeaderFitsReader based on a filename.     
     * @param filename filename
     * @throws URISyntaxException URI syntax problem
     * @throws IOException File not found
     */
    public HeaderFitsReader(final String filename) throws URISyntaxException, IOException {
        final Path path = Paths.get(new URI(filename));
        this.source = Files.newBufferedReader(path, StandardCharsets.UTF_8);
    }

    /**
     * Reads keywords from a reader.
     * @return the list of (keyword,value)
     */
    public List<List<String>> readKeywords() {
        try (BufferedReader reader = new BufferedReader(source)) {
            return reader.lines().filter(line -> line.contains("=")).map(line -> Arrays.asList(line.split(SEPARATOR)))
                         .map(line -> {
                             final String keyword = line.get(0);
                             final String[] valComm = line.get(1).split(" /");
                             final String value = valComm[0].replace("'", "");
                             return Arrays.asList(keyword.trim(), value.trim());
                         }).toList();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
