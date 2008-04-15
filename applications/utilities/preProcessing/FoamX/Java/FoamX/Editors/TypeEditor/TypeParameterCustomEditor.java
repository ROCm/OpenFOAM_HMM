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
package FoamX.Editors.TypeEditor;

import java.awt.Component;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.util.EventObject;

import FoamX.App;
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntryCompoundPanel;

class TypeParameterCustomEditor
    implements TableCellEditor
{
    //--------------------------------------------------------------------------

    protected DictionaryEntryCompoundPanel panel_;
    protected EventListenerList listenerList_;
    protected boolean bAcceptEdit_;
    protected TypeParameter typeParam_;
    transient protected ChangeEvent changeEvent_;
    transient protected ActionEvent actionEvent_;

    //--------------------------------------------------------------------------
    /** TypeParameterCustomEditor constructor. */

    public TypeParameterCustomEditor()
    {
        listenerList_ = new EventListenerList();
        panel_        = new DictionaryEntryCompoundPanel();
        bAcceptEdit_  = false;

        // Subscribe to the edit button's action event so that we can forward
        // the event to the action event listeners of this object.
        panel_.EditButton().addActionListener
        (
            new java.awt.event.ActionListener()
            {
                public void actionPerformed(java.awt.event.ActionEvent evt)
                {
                    // The edit button has been pressed. Inform any interested parties.
                    if (bAcceptEdit_) fireEditActionEvent();
                    bAcceptEdit_ = true;
                }
            }
        );
    }

    //--------------------------------------------------------------------------
    //---- TableCellEditor Interface methods
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
            // Incoming object is a TypeParameter object.
            typeParam_ = (TypeParameter)value;

            // Set the panel text.
            panel_.Label().setText(typeParam_.getValueString());

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
    //---- CellEditor Interface Methods
    //--------------------------------------------------------------------------

    public Object getCellEditorValue()
    {
        return typeParam_;
    }

    //--------------------------------------------------------------------------

    public boolean stopCellEditing()
    {
        typeParam_ = null;
        fireEditingStopped();
        return true;
    }

    //--------------------------------------------------------------------------

    public void cancelCellEditing()
    {
        typeParam_ = null;
        fireEditingCanceled();
    }

    //--------------------------------------------------------------------------

    public boolean isCellEditable(EventObject anEvent)
    {
        // Always editable.
        return true;
    }

    //--------------------------------------------------------------------------

    public boolean shouldSelectCell(EventObject anEvent)
    {
        boolean select = false;

        // Called at the start of the editing process.
        if (anEvent instanceof MouseEvent)
        {
            MouseEvent mevt = (MouseEvent)anEvent;
            if (mevt.getID() == MouseEvent.MOUSE_PRESSED && mevt.getClickCount()>= 2)
            {
                bAcceptEdit_ = false;
                select        = true;
            }
        }

        // Make sure the panel is repainted. Don't know why
        // but sometimes the panel is not rendererd properly
        // despite the table being in edit mode!!
        panel_.repaint();

        return select;
/*
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
*/
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
        for (int i = listeners.length - 2; i>= 0; i -= 2)
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
        // Guaranteed to return a non-null array.
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
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}



