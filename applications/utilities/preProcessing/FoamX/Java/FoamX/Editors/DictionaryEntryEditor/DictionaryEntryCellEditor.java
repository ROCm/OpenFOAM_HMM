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

import java.awt.Component;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.util.EventObject;

import FoamX.App;

public class DictionaryEntryCellEditor
    implements TableCellEditor
{
    //--------------------------------------------------------------------------

    // DictionaryEntry object being edited.
    private DictionaryEntry dictEntry_;
    // Editor being used to edit the current value.
    private TableCellEditor currEditor_;
    // Default cell editor object.
    private DefaultCellEditor defaultEditor_;

    //--------------------------------------------------------------------------
    /** DictionaryEntryCellEditor constructor. */
    public DictionaryEntryCellEditor()
    {
        
        dictEntry_     = null;
        currEditor_    = null;
        defaultEditor_ = null;
    }

    //--------------------------------------------------------------------------
    //---- TableCellEditor Interface Methods
    //--------------------------------------------------------------------------
    /** Implementation of TableCellEditor interface's
      * getTableCellEditorComponent function.
      */
    public Component getTableCellEditorComponent
    (
        JTable  table,
        Object  value,
        boolean isSelected,
        int     row,
        int     column
    )
    {
        Component comp = null;

        try
        {
            // Incoming object must support the DictionaryEntry interface.
            dictEntry_ = (DictionaryEntry)value;

            // See if the DictionaryEntry has a custom editor.
            if (dictEntry_.getEditor() != null)
            {
                currEditor_ = dictEntry_.getEditor();
            }
            else
            {
                // This item does not have a custom editor so use the
                // string-based default.
                if (defaultEditor_ == null)
                {
                    // Create on demand.
                    defaultEditor_ =
                        new DefaultCellEditor(new JTextField());
                }
                currEditor_ = defaultEditor_;
            }

            comp = currEditor_.getTableCellEditorComponent
            (
                table,
                value,
                isSelected,
                row,
                column
            );
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return comp;
    }

    //--------------------------------------------------------------------------
    // CellEditor interface functions.
    //--------------------------------------------------------------------------

    public Object getCellEditorValue()
    {
        // Intercept the "value" object from the editor.
        // Of class String for DefaultCellEditor (JTextField) class.
        // May be null for other editors.

        try
        {
            Object objValue = currEditor_.getCellEditorValue();

            // Update the DictionaryEntry object's value.
            dictEntry_.updateValue(objValue);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }


        // Return the DictionaryEntry object being edited.
        // This is the value that gets sent to the model.
        return dictEntry_;
    }

    //--------------------------------------------------------------------------

    public boolean stopCellEditing()
    {
        if (currEditor_ != null)
        {
            return currEditor_.stopCellEditing();
        }
        return true;
    }

    //--------------------------------------------------------------------------

    public void cancelCellEditing()
    {
        if (currEditor_ != null)
        {
            currEditor_.cancelCellEditing();
        }
    }

    //--------------------------------------------------------------------------

    public boolean isCellEditable(EventObject anEvent)
    {
        // For mouse events, only initiate the editing process is the user
        // has double clicked.
        //if (anEvent instanceof MouseEvent)
        //{
        //    MouseEvent evt = (MouseEvent)anEvent;
        //    return (evt.getID()
        //       == MouseEvent.MOUSE_PRESSED && evt.getClickCount()>= 2);
        //}
        return true;
    }

    //--------------------------------------------------------------------------

    public boolean shouldSelectCell(EventObject anEvent)
    {
        // Called at the start of the editing process.
        if (currEditor_ != null)
        {
            return currEditor_.shouldSelectCell(anEvent);
        }
        return true;
    }

    //--------------------------------------------------------------------------

    public void addCellEditorListener(CellEditorListener l)
    {
        if (currEditor_ != null)
        {
            currEditor_.addCellEditorListener(l);
        }
    }

    //--------------------------------------------------------------------------

    public void removeCellEditorListener(CellEditorListener l)
    {
        if (currEditor_ != null)
        {
            currEditor_.removeCellEditorListener(l);
        }

        // Editing has finished. Reset the references.
        currEditor_ = null;
        dictEntry_  = null;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}




