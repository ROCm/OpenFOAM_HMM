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

import java.awt.event.ActionEvent;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;

import FoamX.App;
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntry;
import FoamX.Editors.DictionaryEntryEditor.EntryCache.DictionaryEntryCache;
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntryCellEditor;
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntryCellRenderer;
import FoamX.Editors.DictionaryEntryEditor.DictionaryNameCellRenderer;
import FoamX.Editors.DictionaryEntryEditor.SelectionEntryCellRenderer;
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntryTableModel;
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntryCompoundEditor;
import FoamX.Editors.TypeEditor.TypeDescriptorCache;

import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;
import FoamXServer.IDictionaryEntry;

public class SelectionEditor extends javax.swing.JDialog
{
    //--------------------------------------------------------------------------

    private static final int DEFAULT_HEIGHT = 400;

    //--------------------------------------------------------------------------

    private DictionaryEntryTableModel tableModel_;
    private IDictionaryEntry dictEntry_;

    //--------------------------------------------------------------------------
    /** Creates new form SelectionEditor */

    public SelectionEditor(java.awt.Frame parent, IDictionaryEntry entry)
    {
        super(parent, "Selection Editor", true);

        dictEntry_ = entry;

        // Initialise the table model.
        tableModel_ = new DictionaryEntryTableModel(entry);

        initComponents();

        // Use 'standard' renderer for name
        table_.getColumn
        (
            tableModel_.getNameColumnTitle()
        ).setCellRenderer(new DictionaryNameCellRenderer());

        // Use specialized renderer for selection entry
        table_.getColumn
        (
            tableModel_.getValueColumnTitle()
        ).setCellRenderer(new SelectionEntryCellRenderer());

        // Use standard dictionaryEntry editor
        table_.getColumn
        (
            tableModel_.getValueColumnTitle()
        ).setCellEditor(new DictionaryEntryCellEditor());


        // Highlight current selection
        int selectedRow = dictEntry_.selection();

        table_.getSelectionModel().setLeadSelectionIndex(selectedRow);
        table_.getSelectionModel().setSelectionInterval
        (
            selectedRow,
            selectedRow
        );
        // Mark selected dictionaryEntry itself
        DictionaryEntryCache subEntry =
            (DictionaryEntryCache)tableModel_.getValueAt
            (
                 selectedRow,
                 DictionaryEntryTableModel.VALUE_COLUMN_INDEX
            );
        subEntry.setCurrent(true);


        // Subscribe to the list selection model event to track changes
        // to current selection
        table_.getSelectionModel().addListSelectionListener
        (
            new ListSelectionListener()
            {
                public void valueChanged(ListSelectionEvent evt)
                {
                    listSelectionChanged(evt);
                }
            }
        );

        // Resize to show the defined rows only.
        // Calculate size: table + panel + extra(=1 row height)
        int height =
            table_.getRowHeight() * (tableModel_.getRowCount() + 1)
            + (int)getMinimumSize().height;  // panel size

        height = Math.min(height, DEFAULT_HEIGHT);

        setSize(new java.awt.Dimension(getSize().width, height));
    }


    //--------------------------------------------------------------------------

    protected void listSelectionChanged(ListSelectionEvent evt)
    {
        // Update the entry status label.
        if (!evt.getValueIsAdjusting())
        {
            int oldRow = dictEntry_.selection();
            int selRow = table_.getSelectionModel().getMinSelectionIndex();

            if (selRow != oldRow)
            {
                // Unset old selected entry.
                DictionaryEntryCache subEntry =
                    (DictionaryEntryCache)tableModel_.getValueAt
                    (
                         oldRow,
                         DictionaryEntryTableModel.VALUE_COLUMN_INDEX
                    );
                subEntry.setCurrent(false);

                // Set the currently selected subEntry.
                subEntry =
                    (DictionaryEntryCache)tableModel_.getValueAt
                    (
                         selRow,
                         DictionaryEntryTableModel.VALUE_COLUMN_INDEX
                    );
                subEntry.setCurrent(true);

                // Mark selected row on this level.
                dictEntry_.selection(selRow);

//                 if (entry != null)
//                 {
//                     TypeDescriptorCache typeDescriptor =
//                         entry.getTypeDescriptor();
//                     System.out.println
//                     (
//                         "even got a description:"
//                         + typeDescriptor.toString()
//                     );
//                     System.out.println
//                     (
//                         "num sub elements"
//                         + entry.getDictEntry().subElements().length
//                     );
//                    TypeDescriptorCache typeDescriptor =
//                        entry.getTypeDescriptor();
//                    if
//                    (
//                        typeDescriptor.isCompound()
//                     && entry.subElements().length == 1)
//                }
            }
        }
    }


    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents()//GEN-BEGIN:initComponents
    {
        buttonPanel_ = new javax.swing.JPanel();
        closeButton_ = new javax.swing.JButton();
        mainPanel_ = new javax.swing.JPanel();
        tableScrollPane_ = new javax.swing.JScrollPane();
        table_ = new javax.swing.JTable();
        
        getContentPane().setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        setName("Shape Editor");
        setModal(true);
        addWindowListener(new java.awt.event.WindowAdapter()
        {
            public void windowClosing(java.awt.event.WindowEvent evt)
            {
                closeDialog(evt);
            }
        });
        
        buttonPanel_.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.RIGHT));
        
        buttonPanel_.setBorder(new javax.swing.border.EtchedBorder());
        closeButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        closeButton_.setText("Close");
        closeButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                closeButtonActionPerformed_(evt);
            }
        });
        
        buttonPanel_.add(closeButton_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints1.weightx = 1.0;
        getContentPane().add(buttonPanel_, gridBagConstraints1);
        
        mainPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        mainPanel_.setBorder(new javax.swing.border.EmptyBorder(new java.awt.Insets(10, 10, 10, 10)));
        table_.setModel(tableModel_);
        table_.setFont(new java.awt.Font("Dialog", 0, 10));
        tableScrollPane_.setViewportView(table_);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints2.weightx = 1.0;
        gridBagConstraints2.weighty = 1.0;
        mainPanel_.add(tableScrollPane_, gridBagConstraints2);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        getContentPane().add(mainPanel_, gridBagConstraints1);
        
        pack();
        java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        setSize(new java.awt.Dimension(350, 300));
        setLocation((screenSize.width-350)/2,(screenSize.height-300)/2);
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

  private void shapeTypeComboActionPerformed (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_shapeTypeComboActionPerformed
  }//GEN-LAST:event_shapeTypeComboActionPerformed

    //--------------------------------------------------------------------------

  private void closeButtonActionPerformed_ (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_closeButtonActionPerformed_
        // Update the dictionary entry and close the dialog.

        // Commit any pending edits.
        table_.editingStopped(new ChangeEvent(this));

        // Close this dialog.
        setVisible(false);
        dispose();
  }//GEN-LAST:event_closeButtonActionPerformed_

    //--------------------------------------------------------------------------
    /** Closes the dialog */
  private void closeDialog(java.awt.event.WindowEvent evt)
    {//GEN-FIRST:event_closeDialog
        setVisible(false);
        dispose();
  }//GEN-LAST:event_closeDialog

    //--------------------------------------------------------------------------

  // Variables declaration - do not modify//GEN-BEGIN:variables
  private javax.swing.JPanel buttonPanel_;
  private javax.swing.JButton closeButton_;
  private javax.swing.JPanel mainPanel_;
  private javax.swing.JScrollPane tableScrollPane_;
  private javax.swing.JTable table_;
  // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}



