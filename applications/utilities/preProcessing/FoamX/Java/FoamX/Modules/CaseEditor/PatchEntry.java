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
package FoamX.Modules.CaseEditor;

import javax.swing.event.*;
import javax.swing.table.*;

import FoamX.App;
import FoamX.Editors.DictionaryEntryEditor.CompoundEntryEvent;
import FoamX.Editors.DictionaryEntryEditor.CompoundEntryListener;
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntry;
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntryCompoundEditor;
import FoamX.Exceptions.FoamXException;

import FoamXServer.FoamXType;
import FoamXServer.IDictionaryEntry;

public class PatchEntry
    implements DictionaryEntry
{
    //--------------------------------------------------------------------------

    protected String patchName_;
    protected IDictionaryEntry dictEntry_;
    protected EventListenerList listenerList_;
    protected DictionaryEntryCompoundEditor compEdit_;

    //--------------------------------------------------------------------------
    /** PatchEntry constructor. */
    public PatchEntry(String patchName, IDictionaryEntry dictEntry)
    {
        try
        {
            // Store IDictionaryEntry reference and patch name.
            dictEntry_ = dictEntry;
            patchName_ = patchName;

            // Make sure that the incoming dictionary entry is a dictionary.
            FoamXType type = dictEntry_.typeDescriptor().type();
            if (type.value() != FoamXType._Type_Dictionary)
            {
                throw new FoamXException("Invalid type for PatchEntry object.");
            }

            // Initialise event listener list.
            listenerList_ = new EventListenerList();

            // Initialise the editor.
            compEdit_ = new DictionaryEntryCompoundEditor();

            // Subscribe to the edit button's action event so that we can
            // invoke the appropriate edit action for the compound types.
            compEdit_.addActionListener
            (
                new java.awt.event.ActionListener()
                {
                    public void actionPerformed(java.awt.event.ActionEvent evt)
                    {
                        fireOpenSubDictionaryEvent();
                    }
                }
            );

        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

    }
    //--------------------------------------------------------------------------

    protected boolean fireOpenSubDictionaryEvent()
    {
        // Send the patch name to the module.

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();

        // Create the event object.
        CompoundEntryEvent event =
            new CompoundEntryEvent
            (
                this,
                CompoundEntryEvent.TYPE_PATCH,
                patchName_
            );

        // Fire the event.
        for (int i = listeners.length - 2; i>= 0; i -= 2)
        {
            if (listeners[i] == CompoundEntryListener.class)
            {
                return
                (
                    (CompoundEntryListener)listeners[i+1]
                ).openSubDictionary(event);
            }
        }

        return false;
    }

    //--------------------------------------------------------------------------
    /**
     * Return an appropriate string representation of the current value.
     * Used by the default cell renderer to set the text for the label.
     */
    public String toString()
    {
        return getEntryValueString();
    }

    //--------------------------------------------------------------------------
    //---- DictionaryEntry Interface Methods
    //--------------------------------------------------------------------------

    public String getEntryName()
    {
        return patchName_;
    }

    public String getEntryDescription()
    {
        return patchName_;
    }

    public String getEntryValueString()
    {
        return "Parameters...";
    }

    public boolean isEditable()
    {
        return dictEntry_.typeDescriptor().editable();
    }

    public boolean isCurrent()
    {
        return false;
    }

    public TableCellRenderer getRenderer()
    {
        return null;
    }

    public TableCellEditor getEditor()
    {
        return compEdit_;
    }

    //--------------------------------------------------------------------------
    /**
     *  Called by the DictionaryCellEditor's getCellEditorValue
     *  method to update the value of this item after editing.
     */
    public boolean updateValue(Object value)
    {
        return false;
    }

    //--------------------------------------------------------------------------

    public void addCompoundEntryListener(CompoundEntryListener l)
    {
        listenerList_.add(CompoundEntryListener.class, l);
    }

    //--------------------------------------------------------------------------

    public void removeCompoundEntryListener(CompoundEntryListener l)
    {
        listenerList_.remove(CompoundEntryListener.class, l);
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}


