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

package FoamX.Editors.DictionaryEntryEditor;

import javax.swing.table.*;
import javax.swing.event.*;
import java.awt.event.*;

import FoamX.App;
import FoamXServer.IDictionaryEntry;
import FoamXServer.ITypeDescriptor;
import FoamXServer.FoamXType;

public class DictionaryEditorPanel
    extends javax.swing.JPanel
{
    //--------------------------------------------------------------------------

    protected static final int  DEFAULT_WIDTH = 400;
    protected DictionaryEntryTableModel tableModel_;
    protected String tableName_;

    //--------------------------------------------------------------------------
    /** Creates a new DictionaryEditorPanel object. */
    public DictionaryEditorPanel(IDictionaryEntry compoundEntry, String key)
    {
        try
        {
            // Initialise the table model.
            tableModel_ = new DictionaryEntryTableModel(compoundEntry);

            IDictionaryEntry[] subEntries = compoundEntry.subElements();

            //tableName_  = compoundEntry.typeDescriptor().name();
            tableName_  = key;
            int type    = compoundEntry.typeDescriptor().type().value();

            tableModel_.setDictionaryKey(key);

            // Initialise the table.
            initComponents();

            // Set single area selection mode.
            parameterTable_.setSelectionMode(0);

            // Use 'standard' renderer for name
            parameterTable_.getColumn
            (
                tableModel_.getNameColumnTitle()
            ).setCellRenderer(new DictionaryNameCellRenderer());

            // Set cell renderer and editor for value column.
            parameterTable_.getColumn
            (
                tableModel_.getValueColumnTitle()
            ).setCellRenderer(new DictionaryEntryCellRenderer());
            parameterTable_.getColumn
            (
                tableModel_.getValueColumnTitle()
            ).setCellEditor(new DictionaryEntryCellEditor());

            // Subscribe to the list selection model event.
            parameterTable_.getSelectionModel().addListSelectionListener
            (
                new ListSelectionListener()
                {
                    public void valueChanged(ListSelectionEvent evt)
                    {
                        listSelectionChanged(evt);
                    }
                }
            );

            // Select the first entry.
            parameterTable_.clearSelection();
            parameterTable_.getSelectionModel().setLeadSelectionIndex(0);

            // Resize to show the defined rows only.
            // Calculate size: table + panel + extra(=1 row height)
            int height =
                parameterTable_.getRowHeight()
                * (tableModel_.getRowCount() + 2)
                + (int)getMinimumSize().getHeight();  // panel size
 
            setPreferredSize(new java.awt.Dimension(DEFAULT_WIDTH, height));
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Provide access to the table model object. */
    public DictionaryEntryTableModel getTableModel()
    {
        return tableModel_;
    }

    //--------------------------------------------------------------------------
    /** Used to set the title of the host internal frame window. */
    public String toString()
    {
        return tableName_;
    }

    //--------------------------------------------------------------------------

    protected void listSelectionChanged(ListSelectionEvent evt)
    {
        // Update the entry status label.
        if (!evt.getValueIsAdjusting())
        {
            // Reset the entry status label.
            statusLabel_.setText("");

            int nSel =
                parameterTable_.getSelectionModel().getMinSelectionIndex();
            if (nSel <tableModel_.getRowCount() && nSel>= 0)
            {
                // Get the currently selected entry.
                DictionaryEntry entry =
                    (DictionaryEntry)tableModel_.getValueAt
                    (
                        nSel,
                        DictionaryEntryTableModel.VALUE_COLUMN_INDEX
                    );
                if (entry != null)
                {
                    statusLabel_.setText(entry.getEntryDescription());
                }
            }
        }
    }

    //--------------------------------------------------------------------------

    public void selectEntry(String entryName)
    {
        try
        {
            // Get index of specified entry.
            int index = tableModel_.getIndexOfEntry(entryName);

            // Give table the focus.
            parameterTable_.requestFocus();

            // Set the specified entry.
            parameterTable_.clearSelection();
            parameterTable_.getSelectionModel().setLeadSelectionIndex(index);
            parameterTable_.getSelectionModel().setSelectionInterval(index, index);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
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
        scrollPane_ = new javax.swing.JScrollPane();
        parameterTable_ = new javax.swing.JTable();
        statusPanel_ = new javax.swing.JPanel();
        statusLabel_ = new javax.swing.JLabel();
        
        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        setBorder(new javax.swing.border.EmptyBorder(new java.awt.Insets(10, 10, 10, 10)));
        setFont(new java.awt.Font("Dialog", 0, 10));
        setAutoscrolls(true);
        scrollPane_.setAutoscrolls(true);
        parameterTable_.setModel(tableModel_);
        parameterTable_.setFont(new java.awt.Font("Dialog", 0, 10));
        parameterTable_.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_ALL_COLUMNS);
        parameterTable_.addKeyListener(new java.awt.event.KeyAdapter()
        {
            public void keyPressed(java.awt.event.KeyEvent evt)
            {
                OnKeyPressed(evt);
            }
        });
        
        scrollPane_.setViewportView(parameterTable_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        add(scrollPane_, gridBagConstraints1);
        
        statusPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        statusPanel_.setBorder(new javax.swing.border.EtchedBorder());
        statusPanel_.setMinimumSize(new java.awt.Dimension(10, 25));
        statusLabel_.setText(" ");
        statusLabel_.setForeground(java.awt.Color.black);
        statusLabel_.setFont(new java.awt.Font("Dialog", 0, 10));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        gridBagConstraints2.weightx = 1.0;
        statusPanel_.add(statusLabel_, gridBagConstraints2);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        add(statusPanel_, gridBagConstraints1);
        
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

  private void OnKeyPressed(java.awt.event.KeyEvent evt)
    {//GEN-FIRST:event_OnKeyPressed
        // Add your handling code here:
        if (evt.getKeyCode() == KeyEvent.VK_TAB)
        {
            int nSel =
                parameterTable_.getSelectionModel().getMinSelectionIndex();
            if (nSel <tableModel_.getRowCount() && nSel>= 0)
            {
                boolean b = parameterTable_.editCellAt(0, 1);
                if (b) App.printMessage("editCellAt true");
            }
        }
  }//GEN-LAST:event_OnKeyPressed

    //--------------------------------------------------------------------------

  // Variables declaration - do not modify//GEN-BEGIN:variables
  private javax.swing.JScrollPane scrollPane_;
  private javax.swing.JTable parameterTable_;
  private javax.swing.JPanel statusPanel_;
  private javax.swing.JLabel statusLabel_;
  // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}



