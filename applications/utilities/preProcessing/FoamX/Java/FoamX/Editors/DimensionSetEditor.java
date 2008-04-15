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

package FoamX.Editors;

import java.awt.*;
import javax.swing.*;
import javax.swing.event.*;
import java.text.*;

import FoamX.App;
import FoamXServer.DimensionSet;
import FoamXServer.DimensionSetHelper;
import FoamXServer.IDictionaryEntry;

public class DimensionSetEditor
    extends javax.swing.JDialog
{
    //--------------------------------------------------------------------------

    protected static final int MASS              = 0;
    protected static final int LENGTH            = 1;
    protected static final int TIME              = 2;
    protected static final int TEMPERATURE       = 3;
    protected static final int MOLES             = 4;
    protected static final int CURRENT           = 5;
    protected static final int LUMINOUSINTENSITY = 6;

    private javax.swing.JLabel[]     labels_     = new javax.swing.JLabel[7];
    private javax.swing.JTextField[] textfields_ = new javax.swing.JTextField[7];
    private IDictionaryEntry entry_;
    private DimensionSet dimensionSet_;

    //--------------------------------------------------------------------------
    /** Constructor for DimensionSetEditor given a non-compound dictionary entry. */

    public DimensionSetEditor(java.awt.Frame parent, IDictionaryEntry entry)
    {
        
        super(parent, "DimensionSet Editor", true);   // Modal.

        // Get the dimension set values from the entry.
        entry_ = entry;
        dimensionSet_ = DimensionSetHelper.extract(entry_.value().value);

        // Common initialisation.
        initialise();

        // Set default dialog position to centre of parent frame.
        if (parent != null)
        {
            setLocation(parent.getLocation().x + (parent.getSize().width  - getSize().width)/2,
                         parent.getLocation().y + (parent.getSize().height - getSize().height)/2);
        }
        else
        {
            Dimension screen = getToolkit().getScreenSize();
            setLocation((screen.getSize().width  - getSize().width)/2,
                         (screen.getSize().height - getSize().height)/2);
        }
    }

    //--------------------------------------------------------------------------

    public DimensionSet getDimensionSet()
    {
        return dimensionSet_;
    }

    //--------------------------------------------------------------------------
    /** Constructor for DimensionSetEditor given a DimensionSet. */

    public DimensionSetEditor(java.awt.Frame parent, DimensionSet dimensionSet)
    {
        super(parent, "DimensionSet Editor", true);

        dimensionSet_ = dimensionSet;

        // Common initialisation.
        initialise();

        // Set default dialog position to centre of parent frame.
        if (parent != null)
        {
            setLocation(parent.getLocation().x + (parent.getSize().width  - getSize().width)/2,
                         parent.getLocation().y + (parent.getSize().height - getSize().height)/2);
        }
        else
        {
            Dimension screen = getToolkit().getScreenSize();
            setLocation((screen.getSize().width  - getSize().width)/2,
                         (screen.getSize().height - getSize().height)/2);
        }
    }

    //--------------------------------------------------------------------------
    /** Common initialisation. */
    public void initialise()
    {
        initComponents();

        // Construct Format object.
        DecimalFormat format = new DecimalFormat();
        format.setGroupingUsed(false);
        format.setMaximumFractionDigits(6);

        // Initialise the dimension edit controls.
        labels_[0] = new javax.swing.JLabel("Mass");
        labels_[1] = new javax.swing.JLabel("Length");
        labels_[2] = new javax.swing.JLabel("Time");
        labels_[3] = new javax.swing.JLabel("Temperature");
        labels_[4] = new javax.swing.JLabel("Moles");
        labels_[5] = new javax.swing.JLabel("Current");
        labels_[6] = new javax.swing.JLabel("Luminous Intensity");

        textfields_[MASS] = new javax.swing.JTextField();
        textfields_[MASS].setText(format.format(dimensionSet_.mass));
        textfields_[LENGTH] = new javax.swing.JTextField();
        textfields_[LENGTH].setText(format.format(dimensionSet_.length));
        textfields_[TIME] = new javax.swing.JTextField();
        textfields_[TIME].setText(format.format(dimensionSet_.time));
        textfields_[TEMPERATURE] = new javax.swing.JTextField();
        textfields_[TEMPERATURE].setText(format.format(dimensionSet_.temperature));
        textfields_[MOLES] = new javax.swing.JTextField();
        textfields_[MOLES].setText(format.format(dimensionSet_.moles));
        textfields_[CURRENT] = new javax.swing.JTextField();
        textfields_[CURRENT].setText(format.format(dimensionSet_.current));
        textfields_[LUMINOUSINTENSITY] = new javax.swing.JTextField();
        textfields_[LUMINOUSINTENSITY].setText(format.format(dimensionSet_.luminousIntensity));

        for (int i=0; i <7; i++)
        {
            labels_[i].setFont(getFont());
            labels_[i].setForeground(java.awt.Color.black);
            textfields_[i].setFont(getFont());
            dimensionPanel_.add(labels_[i]);
            dimensionPanel_.add(textfields_[i]);
        }

        // Resize dialog.
        //pack();
    }

    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents()//GEN-BEGIN:initComponents
    {
        buttonPanel_ = new javax.swing.JPanel();
        btnOK = new javax.swing.JButton();
        btnCancel = new javax.swing.JButton();
        dimensionPanel_ = new javax.swing.JPanel();
        
        setName("Dimension Set Editor");
        setModal(true);
        setFont(new java.awt.Font("Dialog", 0, 10));
        addWindowListener(new java.awt.event.WindowAdapter()
        {
            public void windowClosing(java.awt.event.WindowEvent evt)
            {
                closeDialog(evt);
            }
        });
        
        buttonPanel_.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.RIGHT, 10, 10));
        
        buttonPanel_.setBorder(new javax.swing.border.EtchedBorder());
        buttonPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        btnOK.setFont(new java.awt.Font("Dialog", 0, 10));
        btnOK.setText("OK");
        btnOK.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                btnOKActionPerformed(evt);
            }
        });
        
        buttonPanel_.add(btnOK);
        
        btnCancel.setFont(new java.awt.Font("Dialog", 0, 10));
        btnCancel.setText("Cancel");
        btnCancel.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                btnCancelActionPerformed(evt);
            }
        });
        
        buttonPanel_.add(btnCancel);
        
        getContentPane().add(buttonPanel_, java.awt.BorderLayout.SOUTH);
        
        dimensionPanel_.setLayout(new java.awt.GridLayout(7, 2, 2, 2));
        
        dimensionPanel_.setBorder(new javax.swing.border.EmptyBorder(new java.awt.Insets(10, 10, 10, 10)));
        dimensionPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        getContentPane().add(dimensionPanel_, java.awt.BorderLayout.CENTER);
        
        pack();
        java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        setSize(new java.awt.Dimension(300, 250));
        setLocation((screenSize.width-300)/2,(screenSize.height-250)/2);
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------
    /**
     */
  private void btnCancelActionPerformed (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_btnCancelActionPerformed
        // Shutdown the dialog and return.
        setVisible(false);
        dispose();
    }//GEN-LAST:event_btnCancelActionPerformed

    //--------------------------------------------------------------------------
    private double toDouble(int dimi)
    {
        try
        {
            return Double.parseDouble(textfields_[dimi].getText());
        }
        catch (NumberFormatException ex)
        {
            throw new IllegalArgumentException
            (
                textfields_[dimi].getText()
              + " is not a number for dimension "
              + labels_[dimi].getText()
            );
        }
    }

    //--------------------------------------------------------------------------
    /**
     */
    private void btnOKActionPerformed (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_btnOKActionPerformed

        try
        {
            // Get the values from the text controls and update the dimension set object.
            dimensionSet_.mass              = toDouble(MASS);
            dimensionSet_.length            = toDouble(LENGTH);
            dimensionSet_.time              = toDouble(TIME);
            dimensionSet_.temperature       = toDouble(TEMPERATURE);
            dimensionSet_.current           = toDouble(CURRENT);
            dimensionSet_.moles             = toDouble(MOLES);
            dimensionSet_.luminousIntensity = toDouble(LUMINOUSINTENSITY);
            
            if (entry_ != null)
            {
                // Insert the value into the Any object:
                FoamXServer.FoamXAny value = entry_.value();
                DimensionSetHelper.insert(value.value, dimensionSet_);

                // Update the dictionary entry object.
                entry_.value(value);
            }

            // Close this dialog down.
            setVisible(false);
            dispose();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }//GEN-LAST:event_btnOKActionPerformed

    //--------------------------------------------------------------------------
    /** Closes the dialog */
  private void closeDialog(java.awt.event.WindowEvent evt)
    {//GEN-FIRST:event_closeDialog
        // Shutdown the dialog and return.
        setVisible(false);
        dispose();
  }//GEN-LAST:event_closeDialog

    //--------------------------------------------------------------------------

  // Variables declaration - do not modify//GEN-BEGIN:variables
  private javax.swing.JPanel buttonPanel_;
  private javax.swing.JButton btnOK;
  private javax.swing.JButton btnCancel;
  private javax.swing.JPanel dimensionPanel_;
  // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}


