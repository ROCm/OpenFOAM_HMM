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
package FoamX.Editors.ApplicationEditor;

import java.awt.event.ActionEvent;
import javax.swing.*;
import javax.swing.table.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;

import FoamX.App;

import FoamXServer.CaseServer.IGeometryDescriptor;
import FoamXServer.ITypeDescriptor;

public class FieldsTab
    extends javax.swing.JPanel
{
    //--------------------------------------------------------------------------

    private static final int DEFAULT_WIDTH = 400;

    private ApplicationModel appModel_;
    private DefaultCellEditor stringEditor_;
    private DefaultCellEditor fieldTypeEditor_;
    private DefaultCellEditor geometryTypeEditor_;
    private FieldDimensionCellEditor dimensionEditor_;
    private FieldDimensionCellRenderer dimensionRenderer_;
    private AddFieldAction addFieldAction_;
    private DeleteFieldAction deleteFieldAction_;

    //--------------------------------------------------------------------------
    /** Constructor for PatchPhysicalTypePanel. */
    public FieldsTab(ApplicationModel appModel)
    {
        try
        {
            appModel_ = appModel;

            // Initialise the gui.
            initComponents();

            // Resize table to show the defined rows only.
            int height = (table_.getRowHeight() + table_.getRowMargin()) *
                         (appModel_.getFieldDefinitionModel().getRowCount() + 1);    // Add an extra row to account for the header control.
            table_.setPreferredScrollableViewportSize(new java.awt.Dimension(DEFAULT_WIDTH, height));

            // Initialise the field type, geometry type and dimension editors.
            JComboBox fieldTypeCombo = new JComboBox();
            fieldTypeCombo.setFont(table_.getFont());         // Use same font as table.
            String[] foamTypes = appModel_.getFoamTypes();
            for (int i = 0; i <foamTypes.length; i++)
            {
                ITypeDescriptor fieldType = appModel_.getFoamType(foamTypes[i]);
                fieldTypeCombo.addItem(fieldType.displayName());
            }

            JComboBox geomTypeCombo = new JComboBox();
            geomTypeCombo.setFont(table_.getFont());         // Use same font as table.
            String[] geomTypes = appModel_.getGeometryTypes();
            for (int i = 0; i <geomTypes.length; i++)
            {
                IGeometryDescriptor geomType = appModel_.getGeometryType(geomTypes[i]);
                geomTypeCombo.addItem(geomType.displayName());
            }

            stringEditor_       = new DefaultCellEditor(new JTextField());
            fieldTypeEditor_    = new DefaultCellEditor(fieldTypeCombo);
            geometryTypeEditor_ = new DefaultCellEditor(geomTypeCombo);
            dimensionEditor_    = new FieldDimensionCellEditor();

            // Set cell editor for each column.
            table_.getColumn(FieldDefinitionModel.FIELDNAME_COLUMN).setCellEditor(stringEditor_);
            table_.getColumn(FieldDefinitionModel.FIELDDESC_COLUMN).setCellEditor(stringEditor_);
            table_.getColumn(FieldDefinitionModel.FIELDTYPE_COLUMN).setCellEditor(fieldTypeEditor_);
            table_.getColumn(FieldDefinitionModel.GEOMTYPE_COLUMN).setCellEditor(geometryTypeEditor_);
            table_.getColumn(FieldDefinitionModel.DIMENSION_COLUMN).setCellEditor(dimensionEditor_);

            // Set cell renderer for dimension column.
            dimensionRenderer_ = new FieldDimensionCellRenderer();
            table_.getColumn(FieldDefinitionModel.DIMENSION_COLUMN).setCellRenderer(dimensionRenderer_);

            // Setup the selection listener.
            table_.getSelectionModel().addListSelectionListener
            (
                new ListSelectionListener()
                {
                    public void valueChanged(ListSelectionEvent evt)
                    {
                        if (!table_.getSelectionModel().isSelectionEmpty())
                        {
                            deleteFieldAction_.setEnabled(true);
                        } else {
                            deleteFieldAction_.setEnabled(false);
                        }
                    }
                }
            );

            // Initialise actions and connect to buttons.
            addFieldAction_    = new AddFieldAction();
            deleteFieldAction_ = new DeleteFieldAction();
            addButton_.addActionListener(addFieldAction_);
            //? addFieldAction_.addPropertyChangeListener(new ButtonActionChangedListener(addButton_));
            removeButton_.addActionListener(deleteFieldAction_);
            //? deleteFieldAction_.addPropertyChangeListener(new ButtonActionChangedListener(removeButton_));
            deleteFieldAction_.setEnabled(false);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void updateModel()
    {
        // No update here - App class model extracts info from field model.
    }

    //--------------------------------------------------------------------------
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        subPanel_ = new javax.swing.JPanel();
        scrollPane_ = new javax.swing.JScrollPane();
        table_ = new javax.swing.JTable();
        buttonPanel_ = new javax.swing.JPanel();
        addButton_ = new javax.swing.JButton();
        removeButton_ = new javax.swing.JButton();
        
        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        subPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        scrollPane_.setAutoscrolls(true);
        table_.setModel( appModel_.getFieldDefinitionModel() );
        table_.setFont(new java.awt.Font("Dialog", 0, 10));
        scrollPane_.setViewportView(table_);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTHWEST;
        gridBagConstraints2.weightx = 1.0;
        gridBagConstraints2.weighty = 1.0;
        subPanel_.add(scrollPane_, gridBagConstraints2);
        
        buttonPanel_.setLayout(new javax.swing.BoxLayout(buttonPanel_, javax.swing.BoxLayout.Y_AXIS));
        
        addButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        addButton_.setText("Add");
        addButton_.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        buttonPanel_.add(addButton_);
        
        removeButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        removeButton_.setText("Remove");
        removeButton_.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        buttonPanel_.add(removeButton_);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTHEAST;
        subPanel_.add(buttonPanel_, gridBagConstraints2);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.insets = new java.awt.Insets(10, 10, 10, 10);
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.NORTHWEST;
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        add(subPanel_, gridBagConstraints1);
        
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel subPanel_;
    private javax.swing.JScrollPane scrollPane_;
    private javax.swing.JTable table_;
    private javax.swing.JPanel buttonPanel_;
    private javax.swing.JButton addButton_;
    private javax.swing.JButton removeButton_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------

    class AddFieldAction
    extends AbstractAction
    {
        AddFieldAction()
        {
            putValue(Action.SMALL_ICON, App.getResources().getIcon("AddFieldImage"));
            putValue(Action.NAME, "Add Field");
            putValue(Action.SHORT_DESCRIPTION, "Add Field");
            putValue(Action.LONG_DESCRIPTION, "Add Field");
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                // Add a new field to the table. Easy-Peasy-Lemon-Squeezy.
                appModel_.getFieldDefinitionModel().addNewField();
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    class DeleteFieldAction
    extends AbstractAction
    {
        DeleteFieldAction()
        {
            putValue(Action.NAME, "Delete Field");
            putValue(Action.SMALL_ICON, App.getResources().getIcon("DeleteFieldImage"));
            putValue(Action.SHORT_DESCRIPTION, "Delete Field");
            putValue(Action.LONG_DESCRIPTION, "Delete Field");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Delete the selected dictionary entry.
            try
            {
                // Remove the selected field.
                if (!table_.getSelectionModel().isSelectionEmpty())
                {
                    int minIndex = table_.getSelectionModel().getMinSelectionIndex();
                    int maxIndex = table_.getSelectionModel().getMaxSelectionIndex();
                    // Need to remove from the top so that it doesn't get its knickers in a twist.
                    for (int i=maxIndex; i>= minIndex; i--)
                    {
                        appModel_.getFieldDefinitionModel().removeField(i);
                    }
                    table_.getSelectionModel().clearSelection();
                    deleteFieldAction_.setEnabled(false);
                }
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}


