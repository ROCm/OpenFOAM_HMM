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

import javax.swing.table.*;
import javax.swing.*;
import javax.swing.event.*;

import FoamX.App;
import FoamXServer.CaseServer.IPatchDescriptor;
import FoamXServer.ITypeDescriptor;
//import FoamXServer.CaseServer.IPatchFieldDescriptor;

public class BoundaryDefinitionPanel
    extends javax.swing.JPanel
{
    //--------------------------------------------------------------------------

    private static final int DEFAULT_WIDTH = 100;

    private ApplicationModel appModel_;
    private BoundaryDefinitionModel tableModel_;
    private DefaultCellEditor patchFieldTypeEditor_;
    private JComboBox patchFieldTypeCombo_;

    //--------------------------------------------------------------------------
    /** Constructor for BoundaryDefinitionPanel. */
    public BoundaryDefinitionPanel(ApplicationModel appModel)
    {
        try
        {
            appModel_ = appModel;
            tableModel_    = appModel_.getBoundaryDefinitionModel();

            // Initialise the gui.
            initComponents();

            // Resize table to show the defined rows only.
            int height = (table_.getRowHeight() + table_.getRowMargin()) *
                         (tableModel_.getRowCount() + 1);    // Add an extra row to account for the header control.
            table_.setPreferredScrollableViewportSize(new java.awt.Dimension(DEFAULT_WIDTH, height));

            // Initialise the patch field type editor.
            patchFieldTypeCombo_ = new JComboBox();
            String[] patchFieldTypes = appModel_.getPatchFieldTypes();
            for (int i = 0; i <patchFieldTypes.length; i++)
            {
                ITypeDescriptor patchFieldType = appModel_.getPatchFieldType(patchFieldTypes[i]);
                patchFieldTypeCombo_.addItem(patchFieldType.displayName());
            }
            patchFieldTypeCombo_.setFont(table_.getFont());         // Use same font as table.
            patchFieldTypeEditor_ = new DefaultCellEditor(patchFieldTypeCombo_);

            // Set cell editor for patch field type column.
            table_.getColumn(BoundaryDefinitionModel.PATCH_FIELD_TYPE_COLUMN).setCellEditor(patchFieldTypeEditor_);
            
            // Listen out for changes to the boundary definition model selection.
            tableModel_.addTableModelListener
            (
                new TableModelListener()
                {
                    public void tableChanged(TableModelEvent evt)
                    {
                        boundaryModelChanged(evt);
                    }
                }
            );

            // Initialise the edit controls.
            nameEdit_.setText("");
            displayNameEdit_.setText("");
            descriptionEdit_.setText("");
            nameEdit_.setEnabled(false);
            displayNameEdit_.setEnabled(false);
            descriptionEdit_.setEnabled(false);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    private void boundaryModelChanged(TableModelEvent evt)
    {
        // Another boundary definition has been selected causing the table's data to change.
        patchFieldTypeEditor_.cancelCellEditing();


        BoundaryDefinitionModelItem currItem = tableModel_.getCurrentBoundaryModelItem();
        if (currItem != null && currItem.isBoundaryDefinition())
        {
            nameEdit_.setEnabled(true);
            displayNameEdit_.setEnabled(true);
            descriptionEdit_.setEnabled(true);

            // Update the description text fields.
            nameEdit_.setText(currItem.getName());
            displayNameEdit_.setText(currItem.getDisplayName());
            descriptionEdit_.setText(currItem.getDescription());

            // Update the patch field combo with the valid patch field types for this patch type.
            patchFieldTypeCombo_.removeAllItems();
            IPatchDescriptor patchType = appModel_.getPatchType(currItem.getPatchType());
            String[] patchFieldTypes = patchType.patchFieldTypes();
            for (int i = 0; i <patchFieldTypes.length; i++)
            {
                ITypeDescriptor patchFieldType = appModel_.getPatchFieldType(patchFieldTypes[i]);
                patchFieldTypeCombo_.addItem(patchFieldType.displayName());
            }
        }
        else
        {
            nameEdit_.setText("");
            displayNameEdit_.setText("");
            descriptionEdit_.setText("");
            nameEdit_.setEnabled(false);
            displayNameEdit_.setEnabled(false);
            descriptionEdit_.setEnabled(false);
        }
    }

    //--------------------------------------------------------------------------
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        topPanel_ = new javax.swing.JPanel();
        nameLabel_ = new javax.swing.JLabel();
        nameEdit_ = new javax.swing.JTextField();
        displayNameLabel_ = new javax.swing.JLabel();
        displayNameEdit_ = new javax.swing.JTextField();
        descriptionLabel_ = new javax.swing.JLabel();
        descriptionEdit_ = new javax.swing.JTextField();
        scrollPane_ = new javax.swing.JScrollPane();
        table_ = new javax.swing.JTable();
        
        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        topPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        nameLabel_.setText("Name");
        nameLabel_.setForeground(java.awt.Color.black);
        nameLabel_.setFont(new java.awt.Font("Dialog", 0, 10));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTH;
        topPanel_.add(nameLabel_, gridBagConstraints2);
        
        nameEdit_.setFont(new java.awt.Font("SansSerif", 0, 10));
        nameEdit_.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                nameEditFocusLost(evt);
            }
        });
        
        nameEdit_.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyPressed(java.awt.event.KeyEvent evt) {
                nameEditKeyPressed(evt);
            }
        });
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints2.weightx = 1.0;
        topPanel_.add(nameEdit_, gridBagConstraints2);
        
        displayNameLabel_.setText("Display Name");
        displayNameLabel_.setForeground(java.awt.Color.black);
        displayNameLabel_.setFont(new java.awt.Font("Dialog", 0, 10));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 1;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTH;
        topPanel_.add(displayNameLabel_, gridBagConstraints2);
        
        displayNameEdit_.setFont(new java.awt.Font("SansSerif", 0, 10));
        displayNameEdit_.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                displayNameEditFocusLost(evt);
            }
        });
        
        displayNameEdit_.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyPressed(java.awt.event.KeyEvent evt) {
                displayNameEditKeyPressed(evt);
            }
        });
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 1;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints2.weightx = 1.0;
        topPanel_.add(displayNameEdit_, gridBagConstraints2);
        
        descriptionLabel_.setText("Description");
        descriptionLabel_.setForeground(java.awt.Color.black);
        descriptionLabel_.setFont(new java.awt.Font("Dialog", 0, 10));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 2;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTH;
        topPanel_.add(descriptionLabel_, gridBagConstraints2);
        
        descriptionEdit_.setFont(new java.awt.Font("SansSerif", 0, 10));
        descriptionEdit_.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                descriptionEditFocusLost(evt);
            }
        });
        
        descriptionEdit_.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyPressed(java.awt.event.KeyEvent evt) {
                descriptionEditKeyPressed(evt);
            }
        });
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 2;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints2.weightx = 1.0;
        topPanel_.add(descriptionEdit_, gridBagConstraints2);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints1.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints1.weightx = 1.0;
        add(topPanel_, gridBagConstraints1);
        
        table_.setModel(tableModel_);
        table_.setFont(new java.awt.Font("Dialog", 0, 10));
        table_.setPreferredScrollableViewportSize(new java.awt.Dimension(0, 0));
        scrollPane_.setViewportView(table_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        add(scrollPane_, gridBagConstraints1);
        
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

  private void descriptionEditKeyPressed(java.awt.event.KeyEvent evt)
    {//GEN-FIRST:event_descriptionEditKeyPressed
        if (evt.getKeyCode() == evt.VK_ENTER)
        {
            BoundaryDefinitionModelItem currItem = tableModel_.getCurrentBoundaryModelItem();
            if (currItem != null && currItem.isBoundaryDefinition())
            {
                currItem.setDescription(descriptionEdit_.getText());
                tableModel_.fireBoundaryDefinitionDataChanged();
            }
        }
  }//GEN-LAST:event_descriptionEditKeyPressed

    //--------------------------------------------------------------------------

  private void displayNameEditKeyPressed(java.awt.event.KeyEvent evt)
    {//GEN-FIRST:event_displayNameEditKeyPressed
        if (evt.getKeyCode() == evt.VK_ENTER)
        {
            BoundaryDefinitionModelItem currItem = tableModel_.getCurrentBoundaryModelItem();
            if (currItem != null && currItem.isBoundaryDefinition())
            {
                currItem.setDisplayName(displayNameEdit_.getText());
                tableModel_.fireBoundaryDefinitionDataChanged();
            }
        }
  }//GEN-LAST:event_displayNameEditKeyPressed

    //--------------------------------------------------------------------------

  private void nameEditKeyPressed(java.awt.event.KeyEvent evt)
    {//GEN-FIRST:event_nameEditKeyPressed
        if (evt.getKeyCode() == evt.VK_ENTER)
        {
            BoundaryDefinitionModelItem currItem = tableModel_.getCurrentBoundaryModelItem();
            if (currItem != null && currItem.isBoundaryDefinition())
            {
                currItem.setName(nameEdit_.getText());
                tableModel_.fireBoundaryDefinitionDataChanged();
            }
        }
  }//GEN-LAST:event_nameEditKeyPressed

    //--------------------------------------------------------------------------

  private void nameEditFocusLost (java.awt.event.FocusEvent evt)
    {//GEN-FIRST:event_nameEditFocusLost
        // Update the boundary definition type name.
        BoundaryDefinitionModelItem currItem = tableModel_.getCurrentBoundaryModelItem();
        if (currItem != null && currItem.isBoundaryDefinition())
        {
            currItem.setName(nameEdit_.getText());
            tableModel_.fireBoundaryDefinitionDataChanged();
        }
  }//GEN-LAST:event_nameEditFocusLost

    //--------------------------------------------------------------------------

  private void displayNameEditFocusLost (java.awt.event.FocusEvent evt)
    {//GEN-FIRST:event_displayNameEditFocusLost
        // Update the boundary definition type display name.
        BoundaryDefinitionModelItem currItem = tableModel_.getCurrentBoundaryModelItem();
        if (currItem != null && currItem.isBoundaryDefinition())
        {
            currItem.setDisplayName(displayNameEdit_.getText());
            tableModel_.fireBoundaryDefinitionDataChanged();
        }
  }//GEN-LAST:event_displayNameEditFocusLost

    //--------------------------------------------------------------------------

  private void descriptionEditFocusLost (java.awt.event.FocusEvent evt)
    {//GEN-FIRST:event_descriptionEditFocusLost
        // Update the boundary definition type decription.
        BoundaryDefinitionModelItem currItem = tableModel_.getCurrentBoundaryModelItem();
        if (currItem != null && currItem.isBoundaryDefinition())
        {
            currItem.setDescription(descriptionEdit_.getText());
            tableModel_.fireBoundaryDefinitionDataChanged();
        }
  }//GEN-LAST:event_descriptionEditFocusLost

    //--------------------------------------------------------------------------

  // Variables declaration - do not modify//GEN-BEGIN:variables
  private javax.swing.JPanel topPanel_;
  private javax.swing.JLabel nameLabel_;
  private javax.swing.JTextField nameEdit_;
  private javax.swing.JLabel displayNameLabel_;
  private javax.swing.JTextField displayNameEdit_;
  private javax.swing.JLabel descriptionLabel_;
  private javax.swing.JTextField descriptionEdit_;
  private javax.swing.JScrollPane scrollPane_;
  private javax.swing.JTable table_;
  // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}




