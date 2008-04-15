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

import java.awt.Component;
import java.awt.event.*;
import java.awt.AWTEvent;
import java.text.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.util.EventObject;

import FoamX.App;
import FoamX.Util.FoamXAny;
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntryCompoundPanel;
import FoamX.Editors.DimensionSetEditor;

import FoamXServer.DimensionSet;

public class FieldDimensionCellEditor
    implements TableCellEditor
{
    //--------------------------------------------------------------------------

    protected DimensionSet dimensionSet_;      // DimensionSet object being edited.
    protected DictionaryEntryCompoundPanel panel_;
    protected EventListenerList listenerList_;
    transient protected ChangeEvent changeEvent_;
    transient protected ActionEvent actionEvent_;

    //--------------------------------------------------------------------------

    /** FieldDimensionCellEditor constructor. */
    public FieldDimensionCellEditor()
    {
        dimensionSet_ = null;
        panel_        = new DictionaryEntryCompoundPanel();
        listenerList_ = new EventListenerList();

        // Subscribe to the edit button's action event so that we can forward
        // the event to the action event listeners of this object.
        panel_.EditButton().addActionListener
        (
            new java.awt.event.ActionListener()
            {
                public void actionPerformed(java.awt.event.ActionEvent evt)
                {
                    editButtonActionPerformed(evt);
                }
            }
        );
    }

    //--------------------------------------------------------------------------
    // TableCellEditor interface functions.
    //--------------------------------------------------------------------------
    /** Implementation of TableCellEditor interface's getTableCellEditorComponent function. */
    public Component getTableCellEditorComponent
    (
        JTable  table,
        Object  value,
        boolean isSelected,
        int     row,
        int     column
    )
    {
        try
        {
            dimensionSet_ = (DimensionSet)value;

            // Set the panel text.
            setPanelText();

            // Make sure the label is rendered properly.
            if (isSelected)
            {
                panel_.Label().setForeground(table.getSelectionForeground());
                panel_.Label().setBackground(table.getSelectionBackground());
            } else {
                panel_.Label().setForeground(table.getForeground());
                panel_.Label().setBackground(table.getBackground());
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return panel_;
    }

    //--------------------------------------------------------------------------

    private void setPanelText()
    {
        if (dimensionSet_ != null)
        {
            // Set the panel text.
            panel_.Label().setText(FoamXAny.format(dimensionSet_));
        }
    }

    //--------------------------------------------------------------------------
    // CellEditor interface functions.
    //--------------------------------------------------------------------------

    public Object getCellEditorValue()
    {
        return dimensionSet_;
    }

    //--------------------------------------------------------------------------

    public boolean stopCellEditing()
    {
        fireEditingStopped();
        return true;
    }

    //--------------------------------------------------------------------------

    public void cancelCellEditing()
    {
        fireEditingCanceled();
    }

    //--------------------------------------------------------------------------

    public boolean isCellEditable(EventObject anEvent)
    {
        return true;
    }

    //--------------------------------------------------------------------------

    public boolean shouldSelectCell(EventObject anEvent)
    {
        if (panel_ != null && anEvent instanceof MouseEvent &&
            ((MouseEvent)anEvent).getID() == MouseEvent.MOUSE_PRESSED)
        {
            Component dispatchComponent = SwingUtilities.getDeepestComponentAt(panel_, 3, 3);
            MouseEvent e = (MouseEvent)anEvent;
            MouseEvent e2 = new MouseEvent(dispatchComponent, MouseEvent.MOUSE_RELEASED,
                                            e.getWhen() + 100000, e.getModifiers(), 3, 3, e.getClickCount(),
                                            e.isPopupTrigger());
            dispatchComponent.dispatchEvent(e2);
            e2 = new MouseEvent(dispatchComponent, MouseEvent.MOUSE_CLICKED,
            e.getWhen() + 100001, e.getModifiers(), 3, 3, 1,
            e.isPopupTrigger());
            dispatchComponent.dispatchEvent(e2);
        }
        return false;
    }

    //--------------------------------------------------------------------------

    public void addCellEditorListener(CellEditorListener l)
    {
        listenerList_.add(CellEditorListener.class, l);
    }

    //--------------------------------------------------------------------------

    public void removeCellEditorListener(CellEditorListener l)
    {
        listenerList_.remove(CellEditorListener.class, l);
    }

    //--------------------------------------------------------------------------

    protected void fireEditingStopped()
    {
        Object[] listeners = listenerList_.getListenerList();
        // Process the listeners last to first, notifying those that are interested in this event.
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CellEditorListener.class)
            {
                // Lazily create the event:
                if (changeEvent_ == null)
                {
                    changeEvent_ = new ChangeEvent(this);
                }
                ((CellEditorListener)listeners[i+1]).editingStopped(changeEvent_);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected void fireEditingCanceled()
    {
        // Guaranteed to return a non-null array
        Object[] listeners = listenerList_.getListenerList();
        // Process the listeners last to first, notifying those that are interested in this event.
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CellEditorListener.class)
            {
                // Lazily create the event:
                if (changeEvent_ == null)
                {
                    changeEvent_ = new ChangeEvent(this);
                }
                ((CellEditorListener)listeners[i+1]).editingCanceled(changeEvent_);
            }
        }
    }

    //--------------------------------------------------------------------------

    public void addActionListener(ActionListener l)
    {
        listenerList_.add(ActionListener.class, l);
    }

    //--------------------------------------------------------------------------

    public void removeActionListener(ActionListener l)
    {
        listenerList_.remove(ActionListener.class, l);
    }

    //--------------------------------------------------------------------------

    protected void fireEditActionEvent()
    {
        Object[] listeners = listenerList_.getListenerList();
        // Process the listeners last to first, notifying those that are interested in this event.
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == ActionListener.class)
            {
                // Lazily create the event:
                if (actionEvent_ == null)
                {
                    actionEvent_ = new ActionEvent(this, ActionEvent.ACTION_PERFORMED, "");
                }
                // Fire the event.
                ((ActionListener)listeners[i+1]).actionPerformed(actionEvent_);
            }
        }
    }

    //--------------------------------------------------------------------------

    private void editButtonActionPerformed(java.awt.event.ActionEvent evt)
    {
        // Show the DimensionSet editor dialog.
        DimensionSetEditor dlg = new DimensionSetEditor(null, dimensionSet_);
        dlg.show();
        // Update the panel text.
        setPanelText();
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
