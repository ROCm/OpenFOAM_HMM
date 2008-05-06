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
package FoamX.Editors.DictionaryEntryEditor.EntryCache;

import java.text.*;
import java.awt.Frame;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;

import FoamXServer.CaseServer.*;
import FoamXServer.IDictionaryEntry;
import FoamXServer.IDictionaryEntryHolder;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;
import FoamXServer.FoamXType;

import FoamX.App;
import FoamX.Exceptions.FoamXException;
import FoamX.Editors.CompoundEditor;
import FoamX.Editors.SelectionEditor;
import FoamX.Editors.DimensionSetEditor;
import FoamX.Editors.ListEditor;
import FoamX.Editors.FixedListEditor;
import FoamX.Util.FoamXAny;
import FoamX.Util.FoamXTypeEx;
import FoamX.Editors.TypeEditor.TypeDescriptorCache;
import FoamX.Editors.DictionaryEntryEditor.SelectionEntryCellRenderer;
import FoamX.Editors.DictionaryEntryEditor.CompoundEntryEvent;
import FoamX.Editors.DictionaryEntryEditor.CompoundEntryListener;

public class DictionaryCache
    extends CompoundCache
//    implements DictionaryEntry
{
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    /** DictionaryCache constructor. */
    public DictionaryCache
    (
        IDictionaryEntry dictEntry,
        TypeDescriptorCache typeDescriptor
    )
    {
        super(dictEntry, typeDescriptor);
    }

    //--------------------------------------------------------------------------
    /** DictionaryCache constructor. */
    public DictionaryCache(IDictionaryEntry dictEntry)
    {
        // Invoke the other constructor.
        this
        (
            dictEntry,
            new TypeDescriptorCache(dictEntry.typeDescriptor(), false)
        );
    }

    //--------------------------------------------------------------------------
    /** DictionaryCache constructor. */
    public DictionaryCache(IDictionaryEntry dictEntry, String displayName)
    {
        // Invoke the other constructor.
        this(dictEntry);

        // Use the given display name.
        displayName_ = displayName;
    }

    //--------------------------------------------------------------------------

    public String toString()
    {
        // Initialise if required.
        if (displayValue_.length() == 0)
        {
            displayValue_ = "Dictionary...";
        }
        return displayValue_;
    }

    //--------------------------------------------------------------------------

    protected void editButtonActionPerformed(java.awt.event.ActionEvent evt)
    {
        // The edit button has been pressed. The user wants to edit this entry.
        try
        {
            if (!fireOpenSubDictionaryEvent())
            {
                // Show modal compound editor dialog.
                CompoundEditor editor = new CompoundEditor
                (
                    App.getRootFrame(),
                    dictEntry_
                );
                editor.setTitle(displayName_);
                editor.setVisible(true);
                // Update the cached display string.
                displayValue_ = getCompoundDisplayString(3);
                // Signal that editing has stopped.
                editor_.stopCellEditing();
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    protected boolean fireOpenSubDictionaryEvent()
    {
        // Send the sub-dictionary name to the module.

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();

        // Create the event object.
        CompoundEntryEvent event =
            new CompoundEntryEvent
            (
                this,
                CompoundEntryEvent.TYPE_DICTIONARY,
                typeDescriptor_.getName()
            );

        // Fire the event.
        for (int i = listeners.length - 2; i>= 0; i -= 2)
        {
            if (listeners[i] == CompoundEntryListener.class)
            {
                return
                (
                    (CompoundEntryListener)listeners[i+1])
                        .openSubDictionary(event);
            }
        }

        return false;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}





