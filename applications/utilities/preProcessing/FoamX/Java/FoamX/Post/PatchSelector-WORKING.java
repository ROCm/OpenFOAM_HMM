/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/
package FoamX.Post;

import java.awt.Container;
import java.awt.GridLayout;
import java.awt.Dimension;
import java.awt.event.*;
import java.awt.event.ActionEvent;
import java.awt.Graphics;
import java.awt.GraphicsConfiguration;
import java.io.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;

import com.sun.j3d.utils.universe.*;
import javax.media.j3d.*;
import javax.vecmath.*;

import FoamXServer.*;
import FoamXServer.CasePostServer.ICasePostServer;

import FoamX.App;
import FoamX.CaseManagement.CaseStatusListener;
import FoamX.CaseManagement.CaseStatusEvent;

import FoamX.Post.InteractiveNodes.InteractiveNodeList;
import FoamX.Post.InteractiveNodes.ColoredSurfaceByRefNode;
import FoamX.Post.InteractiveNodes.MeshLinesNode;
import FoamX.Post.InteractiveNodes.ColorBar2DNode;
import FoamX.Post.InteractiveNodes.ArrowsNode;
import FoamX.Post.Shapes.ColoredSurfaceByRef;
import FoamX.Post.Util.Colorizer;
import FoamX.Post.Util.ArrayStats;
import FoamX.Post.Util.PrintManager;
import FoamX.Post.Util.Placement;
import FoamX.Post.Panels.*;

public class PostWindow extends javax.swing.JDialog
{
    //--------------------------------------------------------------------------

    public static final boolean KILL_SERVER = true;

    protected java.awt.Frame parent_;
    protected ICasePostServer casePostServer_;
    protected PostApp pp_;
    protected String caseRoot_;
    protected String caseName_;
    protected EventListenerList listenerList_;
    protected int patchNo_;
    protected int cutNo_;
    protected Transform3D look_;


    //--------------------------------------------------------------------------
    /** Creates new form PostWindow */
    public PostWindow
    (
        java.awt.Frame parent,
        ICasePostServer casePostServer,
        PostApp pp,
        String caseRoot,
        String caseName
    )
    {
        super(parent, "Post Processor", false); // Non modal.

        // Create listener list.
        listenerList_ = new EventListenerList();

        parent_ = parent;
        casePostServer_ = casePostServer;
        pp_ = pp;
        caseRoot_ = caseRoot;
        caseName_ = caseName;

        patchNo_ = 0;
        cutNo_ = 0;

        look_ = new Transform3D();

        initComponents();

        updateButton_ActionPerformed(new ActionEvent(this, 0, "dummy"));

        show();
    }

    //--------------------------------------------------------------------------
    /** Close the case and release the case post server object.
     */
    public void closeCase(boolean killServer)
    {
        try
        {
            // Close down the case server and unlock the case.
            if (killServer && (casePostServer_ != null))
            {
                System.out.println("Killing casePostServer");
                casePostServer_.close();
                casePostServer_ = null;
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }


    //--------------------------------------------------------------------------
    private void initComponents()
    {
        mainPanel_ = new javax.swing.JPanel();

        globalPanel_ = new javax.swing.JPanel();
        updateButton_ = new javax.swing.JButton("Read Times");
        timeBox_ = new javax.swing.JComboBox();
        fieldTypeBox_ = new javax.swing.JComboBox();
        fieldNameBox_ = new javax.swing.JComboBox();
        minField_ = new JTextField("0");
        maxField_ = new JTextField("1");
        vecSizeField_ = new JTextField("1");

        patchPanel_ = new javax.swing.JPanel();
        patchBox_ = new javax.swing.JComboBox();
        createPatchButton_ = new javax.swing.JButton("Create Patch");

        cutPanel_ = new javax.swing.JPanel();
        xField_ = new JTextField("0");
        yField_ = new JTextField("0");
        zField_ = new JTextField("0");
        xNField_ = new JTextField("0");
        yNField_ = new JTextField("0");
        zNField_ = new JTextField("1");
        createCutButton_ = new javax.swing.JButton("Create Cut");
        createMeshCutButton_ = new javax.swing.JButton("Create Mesh Cut");

        viewPanel_ = new javax.swing.JPanel();
        xEyeField_ = new JTextField("1");
        yEyeField_ = new JTextField("0");
        zEyeField_ = new JTextField("0");
        xCenterField_ = new JTextField("0");
        yCenterField_ = new JTextField("0");
        zCenterField_ = new JTextField("0");
        xUpField_ = new JTextField("0");
        yUpField_ = new JTextField("0");
        zUpField_ = new JTextField("1");
        setViewButton_ = new javax.swing.JButton("Set View");

        xResField_ = new JTextField("2000");
        yResField_ = new JTextField("2000");
        printFileName_ = new JFileChooser(new File("."));
        printFileName_.setSelectedFile(new File("FoamX.jpg"));
        printButton_ = new javax.swing.JButton("Print");

        closeButton_ = new javax.swing.JButton("Close");

        // Mainpanel
        mainPanel_.setLayout(new GridLayout(4, 1));
        mainPanel_.add(globalPanel_);
        mainPanel_.add(patchPanel_);
        mainPanel_.add(cutPanel_);
        mainPanel_.add(viewPanel_);
        //mainPanel_.add(closeButton_);

        // Globalpanel
        globalPanel_.setLayout(new GridLayout(7, 2));

        globalPanel_.add(new JLabel("Read Times"));
        globalPanel_.add(updateButton_);

        globalPanel_.add(new JLabel("Time"));
        globalPanel_.add(timeBox_);

        globalPanel_.add(new JLabel("Type"));
        fieldTypeBox_.addItem("volScalarField");
        fieldTypeBox_.addItem("volVectorField");
        globalPanel_.add(fieldTypeBox_);

        globalPanel_.add(new JLabel("Name"));
        globalPanel_.add(fieldNameBox_);

        globalPanel_.add(new JLabel("Min"));
        globalPanel_.add(minField_);

        globalPanel_.add(new JLabel("Max"));
        globalPanel_.add(maxField_);

        globalPanel_.add(new JLabel("VecSize"));
        globalPanel_.add(vecSizeField_);


        // Patchpanel
        patchPanel_.setLayout(new GridLayout(2, 2));
        patchPanel_.add(new JLabel("Patch"));
        patchPanel_.add(patchBox_);

        patchPanel_.add(new JLabel("Create Patch"));
        patchPanel_.add(createPatchButton_);

        // Cutpanel
        cutPanel_.setLayout(new GridLayout(3, 4));
        cutPanel_.add(new JLabel("Point"));
        cutPanel_.add(xField_);
        cutPanel_.add(yField_);
        cutPanel_.add(zField_);

        cutPanel_.add(new JLabel("Normal"));
        cutPanel_.add(xNField_);
        cutPanel_.add(yNField_);
        cutPanel_.add(zNField_);

        cutPanel_.add(new JLabel("Create Cut"));
        cutPanel_.add(createCutButton_);
        cutPanel_.add(new JLabel("Create Mesh Cut"));
        cutPanel_.add(createMeshCutButton_);

        // Viewpanel
        viewPanel_.setLayout(new GridLayout(5,4));

        viewPanel_.add(new JLabel("Eye"));
        viewPanel_.add(xEyeField_);
        viewPanel_.add(yEyeField_);
        viewPanel_.add(zEyeField_);

        viewPanel_.add(new JLabel("Center"));
        viewPanel_.add(xCenterField_);
        viewPanel_.add(yCenterField_);
        viewPanel_.add(zCenterField_);

        viewPanel_.add(new JLabel("Up"));
        viewPanel_.add(xUpField_);
        viewPanel_.add(yUpField_);
        viewPanel_.add(zUpField_);

        viewPanel_.add(new JLabel(""));
        viewPanel_.add(new JLabel(""));
        viewPanel_.add(setViewButton_);
        viewPanel_.add(new JLabel(""));

        viewPanel_.add(xResField_);
        viewPanel_.add(yResField_);
        viewPanel_.add(printButton_);
        viewPanel_.add(closeButton_);

        Container c = getContentPane();
        c.add(mainPanel_);

        // Listeners
        updateButton_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    updateButton_ActionPerformed(evt);
                }
            }
        );
        timeBox_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    timeBox_ActionPerformed(evt);
                }
            }
        );
        fieldTypeBox_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    fieldTypeBox_ActionPerformed(evt);
                }
            }
        );
        createPatchButton_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    createPatchButton_ActionPerformed(evt);
                }
            }
        );
        createCutButton_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    createCutButton_ActionPerformed(evt);
                }
            }
        );
        createMeshCutButton_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    createMeshCutButton_ActionPerformed(evt);
                }
            }
        );
        setViewButton_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    setViewButton_ActionPerformed(evt);
                }
            }
        );
        printButton_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    printButton_ActionPerformed(evt);
                }
            }
        );

        // Window closing
        setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
        addWindowListener
        (
            new java.awt.event.WindowAdapter()
            {
                public void windowClosing(java.awt.event.WindowEvent evt)
                {
                    closeButton_ActionPerformed(new ActionEvent(this, 0, "dummy"));
                }
            }
        );
        closeButton_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    closeButton_ActionPerformed(evt);
                }
            }
        );
        
        

        pack();
        java.awt.Dimension screenSize =
            java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        //setSize(new java.awt.Dimension(300, 200));
        setLocation((screenSize.width-480)/2,(screenSize.height-200)/2);
    }//GEN-END:initComponents


    private void closeButton_ActionPerformed(ActionEvent evt)
    {
        fireClosePostCase
        (
            caseRoot_,
            caseName_
        );
    }

    private void updateButton_ActionPerformed(ActionEvent evt)
    {
        String[] timeSteps = casePostServer_.availableTimeSteps();
        timeBox_.removeAllItems();
        for(int i = 0; i < timeSteps.length; i++)
        {
            timeBox_.addItem(timeSteps[i]);
        }
        timeBox_ActionPerformed(new ActionEvent(this, 0, "dummy"));
    }

    private void setFieldNames(String fieldType)
    {
        try
        {
            if (fieldType == null)
            {
                return;
            }

            String[] fieldNames =
                casePostServer_.getFieldNames(fieldType);

            fieldNameBox_.removeAllItems();
            for(int i = 0; i < fieldNames.length; i++)
            {
                fieldNameBox_.addItem(fieldNames[i]);
            }
        }
        catch (FoamXError fxErr)
        {
            System.out.println("fxErr");
        }
        catch (FoamXIOError ioErr)
        {
            System.out.println("ioErr");
        }
    }

    private void timeBox_ActionPerformed(ActionEvent evt)
    {
        try
        {
            if (timeBox_.getSelectedItem() != null)
            {
                float timeValue =
                    Float.valueOf((String)timeBox_.getSelectedItem()).floatValue();
                int timeIndex = timeBox_.getSelectedIndex();

                casePostServer_.setTime(timeValue, timeIndex);

                String[] patchNames = casePostServer_.getPatchNames();

                patchBox_.removeAllItems();
                for(int i = 0; i < patchNames.length; i++)
                {
                    patchBox_.addItem(patchNames[i]);
                }

                String fieldType = (String)fieldTypeBox_.getSelectedItem();
                setFieldNames(fieldType);
            }
        }
        catch (FoamXError fxErr)
        {
            System.out.println("fxErr");
        }
        catch (FoamXIOError ioErr)
        {
            System.out.println("ioErr");
        }
    }

    private void fieldTypeBox_ActionPerformed(ActionEvent evt)
    {
        String fieldType = (String)fieldTypeBox_.getSelectedItem();
        setFieldNames(fieldType);
    }

    private void createPatchButton_ActionPerformed(ActionEvent evt)
    {
        try
        {
            String fieldType = (String)fieldTypeBox_.getSelectedItem();
            String fieldName = (String)fieldNameBox_.getSelectedItem();
            String patchName = (String)patchBox_.getSelectedItem();

            if
            (
                (fieldName == null)
             || (patchName == null)
             || (fieldType == null)
             || !fieldType.equals("volScalarField")
            )
            {
                return;
            }

            floatListHolder pointsContainer = new floatListHolder();
            longListHolder facesContainer = new longListHolder();
            floatListHolder valuesContainer = new floatListHolder();

            // Read patch geometry + data
            casePostServer_.getTriPatch
            (
                fieldName,
                patchName,
                pointsContainer,
                facesContainer,
                valuesContainer
            );

            float[] points = pointsContainer.value;
            int[] faces = facesContainer.value;
            float[] values = valuesContainer.value;

            if (faces.length == 0)
            {
                System.out.println("Zero size patcht.");
                return;
            }

            float minValue =
                Float.valueOf((String)minField_.getText()).floatValue();
            float maxValue =
                Float.valueOf((String)maxField_.getText()).floatValue();

            if ((minValue == 0.0f) && (maxValue == 0.0f))
            {
                // Autocolor
                ArrayStats stats = new ArrayStats(values);
                minValue = stats.getMin();
                maxValue = stats.getMax();
            }

            String name =
                "patch_" + patchName + "_" + fieldName
              + Integer.toString(patchNo_++);

            // Create geometry object from patch data
            ColoredSurfaceByRefNode ss0 = new ColoredSurfaceByRefNode
            (
                name,
                false,
                false,
                Math.PI / 4.5,
                points,
                null,
                faces,
                values,
                pp_.getColorizer(),
                minValue,
                maxValue
            );
            ColoredSurfaceByRef ss = ss0.getColoredSurfaceByRef();
            ss.setAppearance
            (
                PostApp.createMaterialAppearance()
            );
            ss0.show();
            pp_.getScene().addChild(ss0.getRoot());
            pp_.getShapes().addElement(ss0);
        }
        catch(NumberFormatException nex)
        {
            System.out.println("NumberFormatException");
        }
        catch (FoamXError fxErr)
        {
            System.out.println("fxErr");
        }
        catch (FoamXIOError ioErr)
        {
            System.out.println("ioErr");
        }
    }

    private void createCutButton_ActionPerformed(ActionEvent evt)
    {
        try
        {
            float x = Float.valueOf(xField_.getText()).floatValue();
            float y = Float.valueOf(yField_.getText()).floatValue();
            float z = Float.valueOf(zField_.getText()).floatValue();

            float xN = Float.valueOf(xNField_.getText()).floatValue();
            float yN = Float.valueOf(yNField_.getText()).floatValue();
            float zN = Float.valueOf(zNField_.getText()).floatValue();

            float vecSize = Float.valueOf(vecSizeField_.getText()).floatValue();

            String fieldType = (String)fieldTypeBox_.getSelectedItem();
            String fieldName = (String)fieldNameBox_.getSelectedItem();

            if ((fieldName == null) || (fieldType == null))
            {
                return;
            }
            if ((Math.abs(xN)+Math.abs(yN)+Math.abs(zN)) == 0.0)
            {
                return;
            }

            // Read patch geometry + data
            float[] basePoint = {x, y, z};
            float[] normal = {xN, yN, zN};

            if (fieldType.equals("volScalarField"))
            {
                cutScalar(fieldName, basePoint, normal);
            }
            else if (fieldType.equals("volVectorField"))
            {
                cutVector(fieldName, basePoint, normal, vecSize);
            }
        }
        catch(NumberFormatException nex)
        {
            System.out.println("NumberFormatException");
        }
    }



    private void cutScalar
    (
            String fieldName,
            float[] basePoint,
            float[] normal
    )
    {
        try
        {
            floatListHolder pointsContainer = new floatListHolder();
            longListHolder facesContainer = new longListHolder();
            floatListHolder valuesContainer = new floatListHolder();

            casePostServer_.cutPlane
            (
                fieldName,
                basePoint,
                normal,
                pointsContainer,
                facesContainer,
                valuesContainer
            );
            float[] points = pointsContainer.value;
            int[] faces = facesContainer.value;
            float[] values = valuesContainer.value;

            if (faces.length == 0)
            {
                System.out.println("Zero size cut.");   
                return;
            }
            String name =
                "cut_" + fieldName + "_" + (String)timeBox_.getSelectedItem()
              + "_" + Integer.toString(cutNo_++);

            float minValue =
                Float.valueOf((String)minField_.getText()).floatValue();
            float maxValue =
                Float.valueOf((String)maxField_.getText()).floatValue();

            if ((minValue == 0.0f) && (maxValue == 0.0f))
            {
                // Autocolor
                ArrayStats stats = new ArrayStats(values);
                minValue = stats.getMin();
                maxValue = stats.getMax();
            }

            // Create geometry object from patch data
            ColoredSurfaceByRefNode ss0 = new ColoredSurfaceByRefNode
            (
                name,
                false,
                false,
                Math.PI / 4.5,
                points,
                null,
                faces,
                values,
                pp_.getColorizer(),
                minValue,
                maxValue
            );
            ColoredSurfaceByRef ss = ss0.getColoredSurfaceByRef();
            ss.setAppearance
            (
                PostApp.createMaterialAppearance()
            );
            ss0.show();
            pp_.getScene().addChild(ss0.getRoot());
            pp_.getShapes().addElement(ss0);

            // Add colorBar
            ColorBar2DNode colorBar = new ColorBar2DNode
            (
                "ColorBar_0",
                pp_.getUniverse().getCanvas(),
                new Placement(new Dimension(10, 10), Placement.TOPLEFT),
                new Dimension(20, 128),
                pp_.getColorizer(),
                minValue,
                maxValue
            );
            colorBar.show();
            TransformGroup viewTrans =
                pp_.getUniverse().getViewingPlatform()
                    .getViewPlatformTransform();
            viewTrans.addChild(colorBar.getRoot());
            pp_.getShapes().addElement(colorBar);
            
        }
        catch(NumberFormatException nex)
        {
            System.out.println("NumberFormatException");
        }
        catch(FoamXError fxErr)
        {
            System.out.println("fxErr:" + fxErr);
        }
        catch(FoamXIOError ioErr)
        {
            System.out.println("ioErr:" + ioErr);
        }
    }
        


    private void cutVector
    (
            String fieldName,
            float[] basePoint,
            float[] normal,
            float vecSize
    )
    {
        try
        {
            floatListHolder pointsContainer = new floatListHolder();
            longListHolder facesContainer = new longListHolder();
            floatListHolder valuesContainer = new floatListHolder();

            casePostServer_.cutPlaneVec
            (
                fieldName,
                basePoint,
                normal,
                pointsContainer,
                facesContainer,
                valuesContainer
            );
            float[] points = pointsContainer.value;
            int[] faces = facesContainer.value;
            float[] values = valuesContainer.value;

            if (faces.length == 0)
            {
                System.out.println("Zero size cut.");   
                return;
            }

            // Convert into vectors on face centers.
            int nVecs = values.length / 3;
            float[] sampleCoords = new float[values.length];

            int vertI = 0;
            for(int i = 0; i < nVecs; i++)
            {
                int a = faces[vertI];
                int b = faces[vertI+1];
                int c = faces[vertI+2];

                // Average vertex coords
                float x = points[3*a] + points[3*b] + points[3*c];
                x /= 3.0f;

                float y = points[3*a+1] + points[3*b+1] + points[3*c+1];
                y /= 3.0f;

                float z = points[3*a+2] + points[3*b+2] + points[3*c+2];
                z /= 3.0f;

                sampleCoords[vertI] = x;
                sampleCoords[vertI+1] = y;
                sampleCoords[vertI+2] = z;

                vertI += 3;
            }

//            System.out.println("SampleCoords:");
//            for(int i = 0; i < sampleCoords.length; i+=3)
//            {
//                System.out.println
//                (
//                    "i:" + i + " " + sampleCoords[i] + " " + sampleCoords[i+1]
//                  + " " + sampleCoords[i+2]
//                );
//            }


            String name =
                "cut_" + fieldName + "_" + (String)timeBox_.getSelectedItem()
              + "_" + Integer.toString(cutNo_++);

            // Create geometry object from patch data
            ArrowsNode arr0 = new ArrowsNode
            (
                name,
                sampleCoords,
                values,
                vecSize,
                null
            );
            arr0.show();
            pp_.getScene().addChild(arr0.getRoot());
            pp_.getShapes().addElement(arr0);
        }
        catch(NumberFormatException nex)
        {
            System.out.println("NumberFormatException");
        }
        catch(FoamXError fxErr)
        {
            System.out.println("fxErr:" + fxErr);
        }
        catch(FoamXIOError ioErr)
        {
            System.out.println("ioErr:" + ioErr);
        }
    }




    private void createMeshCutButton_ActionPerformed(ActionEvent evt)
    {
        try
        {
            float x =
                Float.valueOf(xField_.getText()).floatValue();
            float y =
                Float.valueOf(yField_.getText()).floatValue();
            float z =
                Float.valueOf(zField_.getText()).floatValue();

            float xN =
                Float.valueOf(xNField_.getText()).floatValue();
            float yN =
                Float.valueOf(yNField_.getText()).floatValue();
            float zN =
                Float.valueOf(zNField_.getText()).floatValue();

            if ((Math.abs(xN)+Math.abs(yN)+Math.abs(zN)) == 0.0)
            {
                return;
            }

            floatListHolder pointsContainer = new floatListHolder();
            longListHolder edgesContainer = new longListHolder();

            // Read patch geometry + data
            float[] basePoint = {x, y, z};
            float[] normal = {xN, yN, zN};

            casePostServer_.cutMesh
            (
                basePoint,
                normal,
                pointsContainer,
                edgesContainer
            );

            float[] points = pointsContainer.value;
            int[] edges = edgesContainer.value;

            if (edges.length > 0)
            {
                String name =
                    "meshCut_" + Integer.toString(cutNo_++);

                // Create geometry object
                MeshLinesNode ml0 = new MeshLinesNode
                (
                    name,
                    points,
                    edges
                );
                ml0.show();
                pp_.getScene().addChild(ml0.getRoot());
                pp_.getShapes().addElement(ml0);
            }
            else
            {
                System.out.println("Zero size cut.");
            }
        }
        catch(NumberFormatException nex)
        {
            System.out.println("NumberFormatException");
        }
        catch(FoamXError fxErr)
        {
            System.out.println("fxErr:" + fxErr);
        }
        catch(FoamXIOError ioErr)
        {
            System.out.println("ioErr:" + ioErr);
        }
    }

    private void setViewButton_ActionPerformed(ActionEvent evt)
    {
        try
        {
            float xEye =
                Float.valueOf(xEyeField_.getText()).floatValue();
            float yEye =
                Float.valueOf(yEyeField_.getText()).floatValue();
            float zEye =
                Float.valueOf(zEyeField_.getText()).floatValue();

            float xCenter =
                Float.valueOf(xCenterField_.getText()).floatValue();
            float yCenter =
                Float.valueOf(yCenterField_.getText()).floatValue();
            float zCenter =
                Float.valueOf(zCenterField_.getText()).floatValue();

            float xUp =
                Float.valueOf(xUpField_.getText()).floatValue();
            float yUp =
                Float.valueOf(yUpField_.getText()).floatValue();
            float zUp =
                Float.valueOf(zUpField_.getText()).floatValue();

            Point3d eye = new Point3d(xEye, yEye, zEye);
            Point3d center = new Point3d(xCenter, yCenter, zCenter);
            Vector3d up = new Vector3d(xUp, yUp, zUp);

            look_.lookAt(eye, center, up);
            look_.invert();

            ViewingPlatform viewingPlatform =
                pp_.getUniverse().getViewingPlatform();
            TransformGroup viewTrans =
            viewingPlatform.getViewPlatformTransform();
            viewTrans.setTransform(look_);
        }
        catch(NumberFormatException nex)
        {
            System.out.println("NumberFormatException");
        }
    }

    private void printButton_ActionPerformed(ActionEvent evt)
    {
        try
        {
            int xRes =
                Integer.parseInt(xResField_.getText());
            int yRes =
                Integer.parseInt(yResField_.getText());

            System.out.println("xRes:" + xRes + "  yRes:" + yRes);

            // Choose file to print to
            if
            (
                printFileName_.showOpenDialog(this)
             != JFileChooser.APPROVE_OPTION
            )
            {
                return;
            }
            String fileName = printFileName_.getSelectedFile().getName();
            System.out.println("Printing to " + fileName);

            PrintManager printManager = pp_.getPrintManager();
            printManager.print(new Dimension(xRes, yRes), fileName);
        }
        catch(NumberFormatException nex)
        {
            System.out.println("NumberFormatException");
        }

    }


    //--------------------------------------------------------------------------

  // Variables declaration - do not modify//GEN-BEGIN:variables
  private javax.swing.JPanel mainPanel_;
  private javax.swing.JComboBox timeBox_;
  private javax.swing.JComboBox fieldTypeBox_;
  private javax.swing.JComboBox fieldNameBox_;
  private javax.swing.JTextField minField_;
  private javax.swing.JTextField maxField_;
  private javax.swing.JTextField vecSizeField_;
  private javax.swing.JComboBox patchBox_;
  private javax.swing.JButton updateButton_;
  private javax.swing.JButton createPatchButton_;
  private javax.swing.JButton createMeshCutButton_;

  private javax.swing.JPanel globalPanel_;
  private javax.swing.JPanel patchPanel_;
  private javax.swing.JPanel cutPanel_;
  private javax.swing.JTextField xField_;
  private javax.swing.JTextField yField_;
  private javax.swing.JTextField zField_;
  private javax.swing.JTextField xNField_;
  private javax.swing.JTextField yNField_;
  private javax.swing.JTextField zNField_;
  private javax.swing.JButton createCutButton_;
  private javax.swing.JButton closeButton_;

  private javax.swing.JPanel viewPanel_;
  private javax.swing.JTextField xEyeField_;
  private javax.swing.JTextField yEyeField_;
  private javax.swing.JTextField zEyeField_;
  private javax.swing.JTextField xCenterField_;
  private javax.swing.JTextField yCenterField_;
  private javax.swing.JTextField zCenterField_;
  private javax.swing.JTextField xUpField_;
  private javax.swing.JTextField yUpField_;
  private javax.swing.JTextField zUpField_;
  private javax.swing.JButton setViewButton_;

  private javax.swing.JTextField xResField_;
  private javax.swing.JTextField yResField_;
  private javax.swing.JButton printButton_;
  private javax.swing.JFileChooser printFileName_;

  // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    //---- CaseStatusListener Methods
    //--------------------------------------------------------------------------

    public void addCaseStatusListener(CaseStatusListener l)
    {
        listenerList_.add(CaseStatusListener.class, l);
    }

    //--------------------------------------------------------------------------

    public void removeCaseStatusListener(CaseStatusListener l)
    {
        listenerList_.remove(CaseStatusListener.class, l);
    }

    //--------------------------------------------------------------------------

    protected void fireClosePostCase(String caseRoot, String caseName)
    {
        System.out.println("PostWindow: firing closing event for " + caseName);

        // Create event object.
        CaseStatusEvent evt = new CaseStatusEvent(this, caseRoot, caseName);

        // Process the listeners last to first, notifying those that
        // are interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseStatusListener.class)
            {
                ((CaseStatusListener)listeners[i+1]).casePostClosed(evt);
            }
        }
    }


}



