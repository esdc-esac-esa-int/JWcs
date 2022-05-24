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

import io.github.malapert.jwcs.*;
import io.github.malapert.jwcs.proj.AbstractProjection.ProjectionParameter;
import io.github.malapert.jwcs.proj.exception.JWcsException;
import io.github.malapert.jwcs.proj.exception.ProjectionException;
import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JFrame;
import javax.swing.JSlider;

/**
 * ProjectionSelectionPanel lets the user select a projection, applies the
 * selected projection to a group of lines, and displays basic information about
 * the projection.
 *
 * @author Bernhard Jenny, Institute of Cartography, ETH Zurich.
 */
public class ProjectionSelectionPanel extends javax.swing.JPanel {

    /**
     * Logger.
     */
    protected final static Logger LOG = Logger.getLogger(ProjectionSelectionPanel.class.getName());
    private final static long serialVersionUID = 5367008551020527277L;

    /**
     * Countries border.
     */
    private final static String CONTINENTS_PATH = "/io/github/malapert/jwcs/proj/gui/continents.ung";

    /**
     * The lines that are displayed. Must be in geographic coordinates
     * (degrees).
     */
    private List<MapLine> lines = new ArrayList();
    private final List<MapLine> linesFromClient = new ArrayList();
    private final static int STEP_GRID = 5;
    private String previousNameProjection = "";

    /**
     * Creates new form ProjectionSelectionPanel
     */
    public ProjectionSelectionPanel() {
        initComponents();
        Object[] projNames = new Object[]{"ARC", "AZP", "SIN", "STG", "SZP", "TAN", "TPV", "ZEA", "ZPN",
            "AIT", "CAR", "CEA", "CYP", "MER", "MOL", "PAR", "SFL",
            "COD", "COE", "COO", "COP",
            "BON", "PCO","AIR"
        };
        //Object[] projNames = new Object[]{"AIT"};

        projectionComboBox.setModel(new DefaultComboBoxModel(projNames));
    }

    /**
     * Returns a list of latitude lines.
     *
     * @param wcs wcs object
     * @return list of latitude lines
     */
    protected final List<MapLine> drawLatitudeLines(final AbstractJWcs wcs) {
        List<MapLine> latitudes = new ArrayList<>();
        double[] pos1 = new double[2];
        for (int lat = AbstractJWcs.MIN_LATITUDE; lat <= AbstractJWcs.MAX_LATITUDE; lat += STEP_GRID) {
            pos1[0] = Double.NaN;
            for (int lon = AbstractJWcs.MIN_LONGITUDE; lon <= AbstractJWcs.MAX_LONGITUDE; lon += STEP_GRID) {
                double[] pos2;
                try {
                    if (wcs.inside(lon, lat)) {
                        pos2 = wcs.wcs2pix(lon, lat);
                        LOG.log(Level.FINE, "(long,lat)=({0},{1}) --> (x,y)=({2},{3})", new Object[]{lon, lat, pos2[0], pos2[1]});
                        if (wcs.isLineToDraw(pos1, pos2)) {
                            LOG.log(Level.FINE, "plot (x1,y1)=({0},{1}) --> (x2,y2)=({2},{3})", new Object[]{pos1[0], pos1[1], pos2[0], pos2[1]});
                            MapLine line = new MapLine();
                            line.addPoint(pos1[0], pos1[1]);
                            line.addPoint(pos2[0], pos2[1]);
                            latitudes.add(line);
                        }
                        System.arraycopy(pos2, 0, pos1, 0, pos2.length);
                    } else {
                        pos1[0] = Double.NaN;
                    }
                } catch (ProjectionException ex) {
                    pos1[0] = Double.NaN;
                }
            }
        }
        return latitudes;
    }

    /**
     * Returns a list of longitude lines.
     *
     * @param wcs wcs object
     * @return list of latitude lines
     */
    protected final List<MapLine> drawLongitudeLines(final AbstractJWcs wcs) {
        List<MapLine> longitudes = new ArrayList<>();
        double[] pos1 = new double[2];
        for (int lon = AbstractJWcs.MIN_LONGITUDE; lon <= AbstractJWcs.MAX_LONGITUDE; lon += STEP_GRID) {
            pos1[0] = Double.NaN;
            for (int lat = AbstractJWcs.MIN_LATITUDE; lat <= AbstractJWcs.MAX_LATITUDE; lat += STEP_GRID) {
                double[] pos2;
                try {
                    if (wcs.inside(lon, lat)) {
                        pos2 = wcs.wcs2pix(lon, lat);
                        LOG.log(Level.FINE, "(long,lat)=({0},{1}) --> (x,y)=({2},{3})", new Object[]{lon, lat, pos2[0], pos2[1]});
                        if (wcs.isLineToDraw(pos1, pos2)) {
                            LOG.log(Level.FINE, "plot (x1,y1)=({0},{1}) --> (x2,y2)=({2},{3})", new Object[]{pos1[0], pos1[1], pos2[0], pos2[1]});
                            MapLine line = new MapLine();
                            line.addPoint(pos1[0], pos1[1]);
                            line.addPoint(pos2[0], pos2[1]);
                            longitudes.add(line);
                        }
                        System.arraycopy(pos2, 0, pos1, 0, pos2.length);
                    } else {
                        pos1[0] = Double.NaN;
                    }
                } catch (ProjectionException ex) {
                    //LOG.log(Level.WARNING, ex.getMessage());
                    pos1[0] = Double.NaN;
                }
            }
        }
        return longitudes;
    }

    private List<MapLine> drawLines(AbstractJWcs wcs) {
        if (this.linesFromClient.isEmpty()) {
            return new ArrayList<>();
        }
        List<MapLine> projectedLines = new ArrayList<>();
        linesFromClient.stream().map((line) -> line.getPoints()).map((pts) -> {
            MapLine projectedLine = new MapLine();
            pts.stream().filter((pt) -> (wcs.inside(pt.x, pt.y))).forEach((pt) -> {
                try {
                    double[] pos = wcs.wcs2pix(pt.x, pt.y);
                    projectedLine.addPoint(new MapPoint(pos[0], pos[1]));
                } catch (ProjectionException ex) {
                    //Logger.getLogger(ProjectionSelectionPanel.class.getName()).log(Level.SEVERE, null, ex);
                }
            });
            return projectedLine;
        }).forEach((projectedLine) -> {
            projectedLines.add(projectedLine);
        });

        return projectedLines;
    }

    private AbstractJWcs init() throws JWcsException {
        String projName = (String) projectionComboBox.getSelectedItem();
        AbstractJWcs wcs = JWcsMap.getProjection(projName);
        this.lines.clear();
        PV21_Slider.setVisible(false);
        PV21_label.setVisible(false);
        PV21_text.setVisible(false);
        PV22_Slider.setVisible(false);
        PV22_label.setVisible(false);
        PV22_text.setVisible(false);
        PV23_Slider.setVisible(false);
        PV23_label.setVisible(false);
        PV23_text.setVisible(false);
        return wcs;
    }

    private void computeGrid(AbstractJWcs wcs) {
        //this.lines.addAll(drawLatitudeLines(wcs));
        //this.lines.addAll(drawLongitudeLines(wcs));
        this.lines.addAll(drawLines(wcs));
    }

    private void project() {

        try {
            this.errorLabel.setText("");
            JWcsMap jwcsMap = (JWcsMap) init();
            String currentName = jwcsMap.getName();
            ProjectionParameter[] params = jwcsMap.getProjectionParameters();
            switch (params.length) {
                case 1:
                    ProjectionParameter p1 = params[0];
                    double min = p1.getValidInterval()[0];
                    min = Double.isInfinite(min) ? -15 : min;
                    double max = p1.getValidInterval()[1];
                    max = Double.isInfinite(max) ? 15 : max;
                    PV21_Slider.setVisible(true);
                    PV21_Slider.setMinimum((int) min);
                    PV21_Slider.setMaximum((int) max);
                    if (!currentName.equals(previousNameProjection)) {
                        PV21_Slider.setValue((int) p1.getDefaultValue());
                    }
                    PV21_text.setText(p1.getName());
                    PV21_label.setVisible(true);
                    PV21_text.setVisible(true);
                    break;
                case 2:
                    p1 = params[0];
                    min = p1.getValidInterval()[0];
                    min = Double.isInfinite(min) ? -15 : min;
                    max = p1.getValidInterval()[1];
                    max = Double.isInfinite(max) ? 15 : max;
                    PV21_Slider.setVisible(true);
                    PV21_Slider.setMinimum((int) min);
                    PV21_Slider.setMaximum((int) max);
                    //PV21_Slider.setValue((int) p1.getDefaultValue());
                    PV21_text.setText(p1.getName());
                    PV21_label.setVisible(true);
                    PV21_text.setVisible(true);

                    ProjectionParameter p2 = params[1];
                    min = p2.getValidInterval()[0];
                    min = Double.isInfinite(min) ? -15 : min;
                    max = p2.getValidInterval()[1];
                    max = Double.isInfinite(max) ? 15 : max;
                    PV22_Slider.setVisible(true);
                    PV22_Slider.setMinimum((int) min);
                    PV22_Slider.setMaximum((int) max);
                    PV22_text.setText(p2.getName());
                    PV22_label.setVisible(true);
                    PV22_text.setVisible(true);
                    if (!currentName.equals(previousNameProjection)) {
                        PV21_Slider.setValue((int) p1.getDefaultValue());
                        PV22_Slider.setValue((int) p2.getDefaultValue());
                    }
                    break;
                case 3:
                    p1 = params[0];
                    min = p1.getValidInterval()[0];
                    min = Double.isInfinite(min) ? -15 : min;
                    max = p1.getValidInterval()[1];
                    max = Double.isInfinite(max) ? 15 : max;
                    PV21_Slider.setVisible(true);
                    PV21_Slider.setMinimum((int) min);
                    PV21_Slider.setMaximum((int) max);
                    PV21_text.setText(p1.getName());
                    PV21_label.setVisible(true);
                    PV21_text.setVisible(true);

                    p2 = params[1];
                    min = p2.getValidInterval()[0];
                    min = Double.isInfinite(min) ? -15 : min;
                    max = p2.getValidInterval()[1];
                    max = Double.isInfinite(max) ? 15 : max;
                    PV22_Slider.setVisible(true);
                    PV22_Slider.setMinimum((int) min);
                    PV22_Slider.setMaximum((int) max);
                    PV22_text.setText(p2.getName());
                    PV22_label.setVisible(true);
                    PV22_text.setVisible(true);

                    ProjectionParameter p3 = params[2];
                    min = p3.getValidInterval()[0];
                    min = Double.isInfinite(min) ? -15 : min;
                    max = p3.getValidInterval()[1];
                    max = Double.isInfinite(max) ? 15 : max;
                    PV23_Slider.setVisible(true);
                    PV23_Slider.setMinimum((int) min);
                    PV23_Slider.setMaximum((int) max);
                    PV23_text.setText(p3.getName());
                    PV23_label.setVisible(true);
                    PV23_text.setVisible(true);
                    if (!currentName.equals(previousNameProjection)) {
                        PV21_Slider.setValue((int) p1.getDefaultValue());
                        PV22_Slider.setValue((int) p2.getDefaultValue());
                        PV23_Slider.setValue((int) p3.getDefaultValue());
                    }
                    break;
                default:
                    break;

            }
            this.previousNameProjection = currentName;
            jwcsMap.getKeywords().put(AbstractJWcs.CRVAL1, String.valueOf(lon0Slider.getValue()));
            jwcsMap.getKeywords().put(AbstractJWcs.CRVAL2, String.valueOf(lat0Slider.getValue()));
            if (PV21_Slider.isVisible()) {
                jwcsMap.getKeywords().put(AbstractJWcs.PV21, String.valueOf(PV21_Slider.getValue()));
            }
            if (PV22_Slider.isVisible()) {
                jwcsMap.getKeywords().put(AbstractJWcs.PV22, String.valueOf(PV22_Slider.getValue()));
            }
            if (PV23_Slider.isVisible()) {
                jwcsMap.getKeywords().put(AbstractJWcs.PV23, String.valueOf(PV23_Slider.getValue()));
            }
            //jwcsMap.getKeywords().put(AbstractJWcs.PV13, String.valueOf(phip_Slider.getValue()));
            //jwcsMap.getKeywords().put(AbstractJWcs.PV14, String.valueOf(thetap_Slider.getValue()));            
            jwcsMap.doInit();
            computeGrid(jwcsMap);
            map.setLines(lines);
            updateProjectionInfo(jwcsMap);
        } catch (JWcsException ex) {
            String msg = ex.getMessage();
            this.errorLabel.setText(msg);
            //String title = "Error";
            //JOptionPane.showMessageDialog(selectionPanel, msg, title, JOptionPane.ERROR_MESSAGE);
            Logger.getLogger(ProjectionSelectionPanel.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Set the lines that are projected and displayed.
     *
     * @param lines The lines to project. Must be in geographic coordinates
     * (degrees).
     */
    public void setLines(List<MapLine> lines) {
        // store the passed lines
        this.lines = lines;

        // pass the new lines to the map that displays the lines.
        map.setLines(lines);

        // reset the graphical user interface to the Geographical projection.
        projectionComboBox.setSelectedIndex(0);
        project();
    }

    /**
     * Adds lines to the current lines.
     *
     * @param lines the lines to add
     */
    public void addLines(List<MapLine> lines) {
        // store the passed lines
        this.linesFromClient.addAll(lines);
    }

    /**
     * Draw.
     */
    public void draw() {
        projectionComboBox.setSelectedIndex(0);
        project();
    }

    /**
     * Write basic infromation about the projection to the graphical user
     * interface.
     *
     * @param wcs The AbstractProjection that provides the information.
     */
    private void updateProjectionInfo(AbstractJWcs wcs) {
        if (wcs == null) {
            descriptionLabel.setText("-");
        } else {
            descriptionLabel.setText(wcs.getDescription());
        }
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        selectionPanel = new javax.swing.JPanel();
        projectionComboBox = new javax.swing.JComboBox();
        previousProjectionButton = new javax.swing.JButton();
        nextProjectionButton = new javax.swing.JButton();
        map = new io.github.malapert.jwcs.proj.gui.MapComponent();
        infoPanel = new javax.swing.JPanel();
        javax.swing.JLabel descriptionLeadLabel = new javax.swing.JLabel();
        descriptionLabel = new javax.swing.JLabel();
        javax.swing.JLabel longitudeLeadLabel = new javax.swing.JLabel();
        lon0Slider = new javax.swing.JSlider();
        lon0Label = new javax.swing.JLabel();
        jLabel1 = new javax.swing.JLabel();
        lat0Slider = new javax.swing.JSlider();
        lat0Label = new javax.swing.JLabel();
        PV21_Slider = new javax.swing.JSlider();
        PV21_text = new javax.swing.JLabel();
        PV22_Slider = new javax.swing.JSlider();
        PV22_text = new javax.swing.JLabel();
        PV21_label = new javax.swing.JLabel();
        PV22_label = new javax.swing.JLabel();
        PV23_Slider = new javax.swing.JSlider();
        PV23_text = new javax.swing.JLabel();
        PV23_label = new javax.swing.JLabel();
        errorLabel = new javax.swing.JLabel();
        phi0 = new javax.swing.JLabel();
        phip_Slider = new javax.swing.JSlider();
        theta0 = new javax.swing.JLabel();
        thetap_Slider = new javax.swing.JSlider();
        phipLabel = new javax.swing.JLabel();
        thetapLabel = new javax.swing.JLabel();

        //setLayout(new org.netbeans.lib.awtextra.AbsoluteLayout());
        setLayout(new FlowLayout());

        selectionPanel.setPreferredSize(new java.awt.Dimension(100, 40));
        selectionPanel.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 5, 10));

        projectionComboBox.setMaximumRowCount(40);
        projectionComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Plate Carrée (Geographic)", "Cylindrical Equal-Area", "Cylindrical Conformal (Mercator)", "Conical Equidistant", "Conical Equal-Area (Albers)", "Conical Conformal (Lambert)", "Azimuthal Equidistant", "Azimuthal Equal-Area (Lambert)", "Azimuthal Conformal (Stereographic)", "Azimuthal Orthographic", "Sinusoidal", "Pseudoconical Equal-Area (Bonne)" }));
        projectionComboBox.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                projectionComboBoxItemStateChanged(evt);
            }
        });
        selectionPanel.add(projectionComboBox);

        previousProjectionButton.setText("<");
        previousProjectionButton.setPreferredSize(new java.awt.Dimension(50, 29));
        previousProjectionButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                previousProjectionButtonActionPerformed(evt);
            }
        });
        selectionPanel.add(previousProjectionButton);

        nextProjectionButton.setText(">");
        nextProjectionButton.setPreferredSize(new java.awt.Dimension(50, 29));
        nextProjectionButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                nextProjectionButtonActionPerformed(evt);
            }
        });
        selectionPanel.add(nextProjectionButton);

        add(selectionPanel); //, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 0, 604, -1));

        map.setPreferredSize(new java.awt.Dimension(400, 300));
        add(map); //, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 50, 604, 210));

        infoPanel.setBorder(javax.swing.BorderFactory.createCompoundBorder(javax.swing.BorderFactory.createTitledBorder(""), javax.swing.BorderFactory.createEmptyBorder(10, 10, 10, 10)));
        infoPanel.setMinimumSize(new java.awt.Dimension(400, 96));
        infoPanel.setPreferredSize(new java.awt.Dimension(500, 200));

        descriptionLeadLabel.setText("Description");

        descriptionLabel.setText("-");
        descriptionLabel.setMaximumSize(new java.awt.Dimension(300, 16));
        descriptionLabel.setMinimumSize(new java.awt.Dimension(300, 16));
        descriptionLabel.setPreferredSize(new java.awt.Dimension(300, 16));

        longitudeLeadLabel.setText("Longitude of Origin");

        lon0Slider.setMaximum(360);
        lon0Slider.setValue(0);
        lon0Slider.setMinimumSize(new java.awt.Dimension(200, 29));
        lon0Slider.setPreferredSize(new java.awt.Dimension(200, 29));
        lon0Slider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                lon0SliderStateChanged(evt);
            }
        });

        lon0Label.setText("0");
        lon0Label.setMaximumSize(new java.awt.Dimension(50, 16));
        lon0Label.setMinimumSize(new java.awt.Dimension(50, 16));
        lon0Label.setPreferredSize(new java.awt.Dimension(50, 16));

        jLabel1.setText("Latitude of Origin");

        lat0Slider.setMaximum(90);
        lat0Slider.setMinimum(-90);
        lat0Slider.setValue(-90);
        lat0Slider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                lat0SliderStateChanged(evt);
            }
        });

        lat0Label.setText("-90");

        PV21_Slider.setMaximum(10);
        PV21_Slider.setMinimum(-10);
        PV21_Slider.setValue(0);
        PV21_Slider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                PV21_SliderStateChanged(evt);
            }
        });

        PV21_text.setText("mu");

        PV22_Slider.setMaximum(360);
        PV22_Slider.setValue(0);
        PV22_Slider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                PV22_SliderStateChanged(evt);
            }
        });

        PV22_text.setText("gamma");

        PV21_label.setText("0");

        PV22_label.setText("0");

        PV23_Slider.setMaximum(90);
        PV23_Slider.setMinimum(-90);
        PV23_Slider.setValue(90);
        PV23_Slider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                PV23_SliderStateChanged(evt);
            }
        });

        PV23_text.setText("phic");

        PV23_label.setText("90");

        errorLabel.setForeground(new java.awt.Color(230, 17, 17));

        phi0.setText("phip");
        phi0.setEnabled(false);

        phip_Slider.setMaximum(360);
        phip_Slider.setValue(0);
        phip_Slider.setEnabled(false);
        phip_Slider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                phip_SliderStateChanged(evt);
            }
        });

        theta0.setText("thetap");
        theta0.setEnabled(false);

        thetap_Slider.setMaximum(90);
        thetap_Slider.setMinimum(-90);
        thetap_Slider.setToolTipText("");
        thetap_Slider.setValue(90);
        thetap_Slider.setEnabled(false);
        thetap_Slider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                thetap_SliderStateChanged(evt);
            }
        });

        phipLabel.setText("0");
        phipLabel.setEnabled(false);

        thetapLabel.setText("90");
        thetapLabel.setEnabled(false);

        javax.swing.GroupLayout infoPanelLayout = new javax.swing.GroupLayout(infoPanel);
        infoPanel.setLayout(infoPanelLayout);
        infoPanelLayout.setHorizontalGroup(
            infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(infoPanelLayout.createSequentialGroup()
                .addGap(49, 49, 49)
                .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, infoPanelLayout.createSequentialGroup()
                        .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addComponent(PV22_text, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addGap(94, 94, 94))
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addComponent(PV23_text, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
                                .addGap(28, 28, 28))
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addComponent(PV21_text, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addGap(123, 123, 123))
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addComponent(jLabel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addGap(24, 24, 24))
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addComponent(longitudeLeadLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addGap(4, 4, 4))
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(phi0)
                                    .addComponent(theta0))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                        .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(phip_Slider, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(lat0Slider, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(PV21_Slider, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(PV22_Slider, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(PV23_Slider, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(lon0Slider, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(thetap_Slider, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(lon0Label, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(lat0Label, javax.swing.GroupLayout.PREFERRED_SIZE, 26, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(PV21_label)
                                .addComponent(PV22_label)
                                .addComponent(PV23_label)
                                .addComponent(phipLabel))
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addComponent(thetapLabel)
                                .addGap(42, 42, 42)))
                        .addGap(86, 86, 86))
                    .addGroup(infoPanelLayout.createSequentialGroup()
                        .addGap(54, 54, 54)
                        .addComponent(descriptionLeadLabel)
                        .addGap(10, 10, 10)
                        .addComponent(descriptionLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, infoPanelLayout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(errorLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 459, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(89, 89, 89))
        );
        infoPanelLayout.setVerticalGroup(
            infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(infoPanelLayout.createSequentialGroup()
                .addComponent(errorLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 18, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(lat0Label)
                    .addGroup(infoPanelLayout.createSequentialGroup()
                        .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(descriptionLeadLabel)
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addGap(1, 1, 1)
                                .addComponent(descriptionLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addGap(10, 10, 10)
                                .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(longitudeLeadLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(lon0Label, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGap(2, 2, 2))
                            .addComponent(lon0Slider, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(lat0Slider, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addComponent(jLabel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addGap(2, 2, 2)))))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(infoPanelLayout.createSequentialGroup()
                        .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addGap(0, 0, Short.MAX_VALUE)
                                .addComponent(PV21_label))
                            .addComponent(PV21_Slider, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE))
                        .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addGap(10, 10, 10)
                                .addComponent(PV22_label))
                            .addGroup(infoPanelLayout.createSequentialGroup()
                                .addGap(2, 2, 2)
                                .addComponent(PV22_Slider, javax.swing.GroupLayout.PREFERRED_SIZE, 32, Short.MAX_VALUE)
                                .addGap(4, 4, 4))))
                    .addGroup(infoPanelLayout.createSequentialGroup()
                        .addComponent(PV21_text, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(PV22_text, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addGap(4, 4, 4)))
                .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(PV23_label)
                    .addGroup(infoPanelLayout.createSequentialGroup()
                        .addGap(4, 4, 4)
                        .addComponent(PV23_text, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addComponent(PV23_Slider, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(infoPanelLayout.createSequentialGroup()
                        .addComponent(phipLabel)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addComponent(phip_Slider, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
                    .addComponent(phi0, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(infoPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(theta0, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(thetap_Slider, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
                    .addComponent(thetapLabel))
                .addGap(23, 23, 23))
        );

        add(infoPanel); //, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 268, 604, 250));
    }// </editor-fold>//GEN-END:initComponents

    private void lon0SliderStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_lon0SliderStateChanged
        JSlider slider = (JSlider) evt.getSource();
        lon0Label.setText(Integer.toString(slider.getValue()));
        //if (!slider.getValueIsAdjusting()) {
        project();
        //}
    }//GEN-LAST:event_lon0SliderStateChanged

    private void projectionComboBoxItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_projectionComboBoxItemStateChanged
        if (evt.getStateChange() == java.awt.event.ItemEvent.SELECTED) {
            project();
        }
    }//GEN-LAST:event_projectionComboBoxItemStateChanged

    private void previousProjectionButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_previousProjectionButtonActionPerformed
        int id = projectionComboBox.getSelectedIndex() - 1;
        if (id >= 0) {
            projectionComboBox.setSelectedIndex(id);
            project();
        }
    }//GEN-LAST:event_previousProjectionButtonActionPerformed

    private void nextProjectionButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_nextProjectionButtonActionPerformed
        int id = projectionComboBox.getSelectedIndex() + 1;
        if (id < projectionComboBox.getItemCount()) {
            projectionComboBox.setSelectedIndex(id);
            project();
        }
    }//GEN-LAST:event_nextProjectionButtonActionPerformed

    private void lat0SliderStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_lat0SliderStateChanged
        JSlider slider = (JSlider) evt.getSource();
        lat0Label.setText(Integer.toString(slider.getValue()));
        //if (!slider.getValueIsAdjusting()) {
        project();
        //}
    }//GEN-LAST:event_lat0SliderStateChanged

    private void PV21_SliderStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_PV21_SliderStateChanged
        JSlider slider = (JSlider) evt.getSource();
        PV21_label.setText(Integer.toString(slider.getValue()));
        project();
    }//GEN-LAST:event_PV21_SliderStateChanged

    private void PV22_SliderStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_PV22_SliderStateChanged
        JSlider slider = (JSlider) evt.getSource();
        PV22_label.setText(Integer.toString(slider.getValue()));
        project();
    }//GEN-LAST:event_PV22_SliderStateChanged

    private void PV23_SliderStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_PV23_SliderStateChanged
        JSlider slider = (JSlider) evt.getSource();
        PV23_label.setText(Integer.toString(slider.getValue()));
        project();
    }//GEN-LAST:event_PV23_SliderStateChanged

    private void phip_SliderStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_phip_SliderStateChanged
        JSlider slider = (JSlider) evt.getSource();
        this.phipLabel.setText(Integer.toString(slider.getValue()));
        project();
    }//GEN-LAST:event_phip_SliderStateChanged

    private void thetap_SliderStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_thetap_SliderStateChanged
        JSlider slider = (JSlider) evt.getSource();
        this.thetapLabel.setText(Integer.toString(slider.getValue()));
        project();
    }//GEN-LAST:event_thetap_SliderStateChanged

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JSlider PV21_Slider;
    private javax.swing.JLabel PV21_label;
    private javax.swing.JLabel PV21_text;
    private javax.swing.JSlider PV22_Slider;
    private javax.swing.JLabel PV22_label;
    private javax.swing.JLabel PV22_text;
    private javax.swing.JSlider PV23_Slider;
    private javax.swing.JLabel PV23_label;
    private javax.swing.JLabel PV23_text;
    private javax.swing.JLabel descriptionLabel;
    private javax.swing.JLabel errorLabel;
    private javax.swing.JPanel infoPanel;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel lat0Label;
    private javax.swing.JSlider lat0Slider;
    private javax.swing.JLabel lon0Label;
    private javax.swing.JSlider lon0Slider;
    private io.github.malapert.jwcs.proj.gui.MapComponent map;
    private javax.swing.JButton nextProjectionButton;
    private javax.swing.JLabel phi0;
    private javax.swing.JLabel phipLabel;
    private javax.swing.JSlider phip_Slider;
    private javax.swing.JButton previousProjectionButton;
    private javax.swing.JComboBox projectionComboBox;
    private javax.swing.JPanel selectionPanel;
    private javax.swing.JLabel theta0;
    private javax.swing.JLabel thetapLabel;
    private javax.swing.JSlider thetap_Slider;
    // End of variables declaration//GEN-END:variables

    /**
     * Create a window, ask the user for lines to display, import the lines, and
     * display them.
     *
     * @throws java.io.IOException Cannot import countries
     */
    public static void createWindow() throws IOException {
        // create a new window
        JFrame mapWindow = new JFrame("JWcs - Projection");
        mapWindow.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        ProjectionSelectionPanel panel = new ProjectionSelectionPanel();
        mapWindow.getContentPane().add(panel, BorderLayout.CENTER);
        mapWindow.pack();
        mapWindow.setLocationRelativeTo(null); // center on screen
        mapWindow.setVisible(true);
        URL url = ProjectionSelectionPanel.class.getResource(CONTINENTS_PATH);
        InputStream stream = url.openStream();
        List<MapLine> lines = UngenerateImporter.importData(stream);
        // pass the lines to the map component
        panel.addLines(lines);
        panel.draw();
    }
    
    /**
     * Main method.
     * @param args arguments 
     */
    public static void main(String[] args) {
        java.awt.EventQueue.invokeLater(() -> {
            try {
                createWindow();
            } catch (IOException ex) {
                LOG.log(Level.SEVERE, null, ex);
            }
        });
    }    
}
