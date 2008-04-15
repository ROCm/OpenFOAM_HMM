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

import java.util.*;
import java.awt.event.ActionEvent;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;

import FoamX.App;

public class StringListEditor
    extends javax.swing.JDialog
{
    //--------------------------------------------------------------------------

    private static final int DEFAULT_HEIGHT = 400;

    //--------------------------------------------------------------------------

    private ListEntryModel tableModel_;
    private AddElementAction addElementAction_;
    private DeleteElementAction deleteElementAction_;
    private Vector values_;
    private Vector labels_;
    private boolean wasCancelled_;

    //--------------------------------------------------------------------------
    /** StringListEditor constructor */
    public StringListEditor(java.awt.Frame parent, java.util.Vector values)
    {
        
        super(parent, "List Editor", true);   // Modal.

        labels_ = null;
        values_ = values;

        // Initialise the table model.
        tableModel_ = new ListEntryModel(values_);

        // Initialise.
        init();
    }

    //--------------------------------------------------------------------------
    /** StringListEditor constructor */
    public StringListEditor
    (
        java.awt.Frame parent,
        String title,
        java.util.Vector labels,
        java.util.Vector values
    )
    {
        super(parent, title, true);   // Modal.

        labels_ = labels;
        values_ = values;

        // Initialise the table model.
        tableModel_ = new ListEntryModel(labels_, values_);

        // Initialise.
        init();

        // Disable add/remove buttons.
        btnInsertElement_.setVisible(false);
        btnDeleteElement_.setVisible(false);
    }

    //--------------------------------------------------------------------------

    private void init()
    {
        initComponents();

        // Resize to show the defined rows only.
        // Calculate size: table + some extra (= 1 row height) or
        //                 space needed for panel itself (so with empty table)
        int height =
            (
                listTable_.getRowHeight() * (tableModel_.getRowCount() + 1)
              + (int)getMinimumSize().getHeight()
            );    // panel size
        height = Math.min(height, DEFAULT_HEIGHT);

        setSize(new java.awt.Dimension(getWidth(), height));

        // Initialise actions and hookup to buttons.
        addElementAction_    = new AddElementAction();
        deleteElementAction_ = new DeleteElementAction();

        btnInsertElement_.addActionListener(addElementAction_);
        //? addElementAction_.addPropertyChangeListener
        //? (new ButtonActionChangedListener(btnInsertElement_));

        btnDeleteElement_.addActionListener(deleteElementAction_);
        //? deleteElementAction_.addPropertyChangeListener
        //? (new ButtonActionChangedListener(btnDeleteElement_));

        // Listen out for selection events and enable the appropriate actions.
        listTable_.getSelectionModel().addListSelectionListener
        (
            new ListSelectionListener()
            {
                public void valueChanged(ListSelectionEvent evt)
                {
                    // Enabled if selection is not empty.
                    deleteElementAction_.setEnabled
                    (
                        !listTable_.getSelectionModel().isSelectionEmpty()
                    );
                }
            }
        );

        wasCancelled_ = true;
        // Set cell renderer and editor for value column.
        //listTable_.getColumn
        //    (
        //        tableModel_.getValueColumnTitle()
        //    ).setCellRenderer(new DictionaryEntryCellRenderer());
        //listTable_.getColumn
        //    (
        //        tableModel_.getValueColumnTitle()
        //    ).setCellEditor(new DictionaryEntryCellEditor());
    }

    //--------------------------------------------------------------------------

    public boolean wasCancelled()
    {
        return wasCancelled_;
    }

    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents()//GEN-BEGIN:initComponents
    {
        mainPanel_ = new javax.swing.JPanel();
        tableScrollPane_ = new javax.swing.JScrollPane();
        listTable_ = new javax.swing.JTable();
        listButtonPanel_ = new java.awt.Panel();
        btnInsertElement_ = new javax.swing.JButton();
        btnDeleteElement_ = new javax.swing.JButton();
        buttonPanel_ = new javax.swing.JPanel();
        cancelButton_ = new javax.swing.JButton();
        closeButton_ = new javax.swing.JButton();
        
        getContentPane().setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        setName("List Editor");
        setModal(true);
        setFont(new java.awt.Font("Dialog", 0, 10));
        addWindowListener(new java.awt.event.WindowAdapter()
        {
            public void windowClosing(java.awt.event.WindowEvent evt)
            {
                closeDialog(evt);
            }
        });
        
        mainPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        mainPanel_.setBorder(new javax.swing.border.EmptyBorder(new java.awt.Insets(10, 10, 10, 10)));
        listTable_.setModel(tableModel_);
        listTable_.setFont(new java.awt.Font("Dialog", 0, 10));
        tableScrollPane_.setViewportView(listTable_);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTHWEST;
        gridBagConstraints2.weightx = 1.0;
        gridBagConstraints2.weighty = 1.0;
        mainPanel_.add(tableScrollPane_, gridBagConstraints2);
        
        listButtonPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints3;
        
        listButtonPanel_.setFont(new java.awt.Font("Dialog", 0, 11));
        listButtonPanel_.setBackground(new java.awt.Color(204, 204, 204));
        listButtonPanel_.setForeground(java.awt.Color.black);
        btnInsertElement_.setFont(new java.awt.Font("Dialog", 0, 10));
        btnInsertElement_.setText("Insert...");
        gridBagConstraints3 = new java.awt.GridBagConstraints();
        gridBagConstraints3.gridx = 0;
        gridBagConstraints3.gridy = 0;
        gridBagConstraints3.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints3.insets = new java.awt.Insets(0, 0, 5, 0);
        listButtonPanel_.add(btnInsertElement_, gridBagConstraints3);
        
        btnDeleteElement_.setFont(new java.awt.Font("Dialog", 0, 10));
        btnDeleteElement_.setText("Delete...");
        gridBagConstraints3 = new java.awt.GridBagConstraints();
        gridBagConstraints3.gridx = 0;
        gridBagConstraints3.gridy = 1;
        gridBagConstraints3.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints3.insets = new java.awt.Insets(0, 0, 5, 0);
        listButtonPanel_.add(btnDeleteElement_, gridBagConstraints3);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTHEAST;
        mainPanel_.add(listButtonPanel_, gridBagConstraints2);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        getContentPane().add(mainPanel_, gridBagConstraints1);
        
        buttonPanel_.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.RIGHT));
        
        buttonPanel_.setBorder(new javax.swing.border.EtchedBorder());
        cancelButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        cancelButton_.setText("Cancel");
        cancelButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnCancel(evt);
            }
        });
        
        buttonPanel_.add(cancelButton_);
        
        closeButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        closeButton_.setText("OK");
        closeButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnClose(evt);
            }
        });
        
        buttonPanel_.add(closeButton_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints1.weightx = 1.0;
        getContentPane().add(buttonPanel_, gridBagConstraints1);
        
        pack();
        java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        setSize(new java.awt.Dimension(350, 300));
        setLocation((screenSize.width-350)/2,(screenSize.height-300)/2);
    }//GEN-END:initComponents

  private void OnCancel(java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnCancel

        // Set cancelled flag.
        wasCancelled_ = true;

        // Close dialog.
        setVisible(false);
        dispose();

  }//GEN-LAST:event_OnCancel

    //--------------------------------------------------------------------------

  private void shapeTypeComboActionPerformed (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_shapeTypeComboActionPerformed
  }//GEN-LAST:event_shapeTypeComboActionPerformed

    //--------------------------------------------------------------------------

  private void OnClose (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnClose

        // Set cancelled flag.
        wasCancelled_ = false;

        // Commit any pending edits.
        listTable_.editingStopped(new ChangeEvent(this));

        // Close dialog.
        setVisible(false);
        dispose();

  }//GEN-LAST:event_OnClose

    //--------------------------------------------------------------------------
    /** Closes the dialog */
  private void closeDialog(java.awt.event.WindowEvent evt)
    {//GEN-FIRST:event_closeDialog

        wasCancelled_ = true;

        setVisible(false);
        dispose();

  }//GEN-LAST:event_closeDialog

    //--------------------------------------------------------------------------

  // Variables declaration - do not modify//GEN-BEGIN:variables
  private javax.swing.JPanel mainPanel_;
  private javax.swing.JScrollPane tableScrollPane_;
  private javax.swing.JTable listTable_;
  private java.awt.Panel listButtonPanel_;
  private javax.swing.JButton btnInsertElement_;
  private javax.swing.JButton btnDeleteElement_;
  private javax.swing.JPanel buttonPanel_;
  private javax.swing.JButton cancelButton_;
  private javax.swing.JButton closeButton_;
  // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    private static class ListEntryModel
    extends AbstractTableModel
    {
        //----------------------------------------------------------------------

        static final int      PARAMETERNAME_INDEX  = 0;
        static final int      PARAMETERVALUE_INDEX = 1;
        static final String   PARAMETERNAME_TEXT   = "Parameter";
        static final String   PARAMETERVALUE_TEXT  = "Value";
        static final String[] s_columnNames        = {
        PARAMETERNAME_TEXT, PARAMETERVALUE_TEXT };

        java.util.Vector labelVector_;
        java.util.Vector valueVector_;

        //----------------------------------------------------------------------
        /** ListEntryModel constructor for list types. */
        ListEntryModel(java.util.Vector valueVector)
        {
            labelVector_ = null;
            valueVector_ = valueVector;
        }

        //----------------------------------------------------------------------
        /** ListEntryModel constructor for list types. */
        ListEntryModel
        (
            java.util.Vector labelVector,
            java.util.Vector valueVector
        )
        {
            labelVector_ = labelVector;
            valueVector_ = valueVector;
        }

        //----------------------------------------------------------------------
        //---- AbstractTableModel methods
        //----------------------------------------------------------------------

        public int getRowCount()
        {
            return valueVector_.size();
        }

        //----------------------------------------------------------------------

        public int getColumnCount()
        {
            return s_columnNames.length;
        }

        //----------------------------------------------------------------------

        public String getColumnName(int columnIndex)
        {
            return s_columnNames[columnIndex];
        }

        //----------------------------------------------------------------------

        public Class getColumnClass(int columnIndex)
        {
            // Everything is of type String.
            return String.class;
        }

        //----------------------------------------------------------------------

        public boolean isCellEditable(int rowIndex, int columnIndex)
        {
            // Check for value column.
            return (columnIndex == PARAMETERVALUE_INDEX);
        }

        //----------------------------------------------------------------------

        public Object getValueAt(int rowIndex, int columnIndex)
        {
            String value = "";

            switch (columnIndex)
            {
                case PARAMETERNAME_INDEX:
                    if (labelVector_ != null)
                    {
                        value = labelVector_.get(rowIndex).toString();
                    }
                    else
                    {
                        // Return row number.
                        value = Integer.toString(rowIndex);
                    }
                    break;
                case PARAMETERVALUE_INDEX:
                    value = valueVector_.get(rowIndex).toString();
                    break;
                default:
                    // Bugger and damnation.
                    break;
            }

            return value;
        }

        //----------------------------------------------------------------------

        public void setValueAt(Object aValue, int rowIndex, int columnIndex)
        {
            switch (columnIndex)
            {
                case PARAMETERVALUE_INDEX:
                    String value = (String)aValue;
                    valueVector_.set(rowIndex, value);
                    break;
                case PARAMETERNAME_INDEX:
                default:
                    // Bugger.
                    break;
            }
        }

        //----------------------------------------------------------------------
        //----------------------------------------------------------------------
        //----------------------------------------------------------------------

        void addElement(String value)
        {
            try
            {
                // Add new elememt to the list.
                valueVector_.addElement(value);
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }

        //----------------------------------------------------------------------

        void removeElement(int index)
        {
            try
            {
                // Check for valid row index.
                if (index <valueVector_.size())
                {
                    // Remove this entry from the model.
                    valueVector_.removeElementAt(index);
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

    class AddElementAction
    extends AbstractAction
    {
        AddElementAction()
        {
            super("AddElement");
            putValue(Action.NAME, "Add Element");
            putValue(Action.SHORT_DESCRIPTION, "Add Element");
            putValue(Action.LONG_DESCRIPTION, "Add Element");
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                tableModel_.addElement("");
                listTable_.updateUI();
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    class DeleteElementAction
    extends AbstractAction
    {
        DeleteElementAction()
        {
            super("DeleteElement");
            putValue(Action.NAME, "Delete Element");
            putValue(Action.SHORT_DESCRIPTION, "Delete Element");
            putValue(Action.LONG_DESCRIPTION, "Delete Element");
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                // Remove the selected elements.
                if (!listTable_.getSelectionModel().isSelectionEmpty())
                {
                    int minIndex =
                        listTable_.getSelectionModel().getMinSelectionIndex();
                    int maxIndex =
                        listTable_.getSelectionModel().getMaxSelectionIndex();
                    // Need to remove from the top so that it doesn't get its
                    // knickers in a twist.
                    for (int i=maxIndex; i>= minIndex; i--)
                    {
                        tableModel_.removeElement(i);
                    }
                    listTable_.getSelectionModel().clearSelection();
                    deleteElementAction_.setEnabled(false);
                }
                listTable_.updateUI();
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }
}



