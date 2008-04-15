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

public class SelectionCache
    extends CompoundCache
//    implements DictionaryEntry
{
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    /** SelectionCache constructor. */
    public SelectionCache
    (
        IDictionaryEntry dictEntry,
        TypeDescriptorCache typeDescriptor
    )
    {
        super(dictEntry, typeDescriptor);
    }

    //--------------------------------------------------------------------------
    /** SelectionCache constructor. */
    public SelectionCache(IDictionaryEntry dictEntry)
    {
        // Invoke the other constructor.
        this
        (
            dictEntry,
            new TypeDescriptorCache(dictEntry.typeDescriptor(), false)
        );
    }

    //--------------------------------------------------------------------------
    /** SelectionCache constructor. */
    public SelectionCache(IDictionaryEntry dictEntry, String displayName)
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
            displayValue_ = getSelectionDisplayString() + " ...";
        }
        return displayValue_;
    }

    //--------------------------------------------------------------------------

    public void toStringRaw(Vector argListVector)
    {
        getSelectionValueString(argListVector);
    }

    //--------------------------------------------------------------------------

    protected void editButtonActionPerformed(java.awt.event.ActionEvent evt)
    {
        // The edit button has been pressed. The user wants to edit this entry.
        try
        {
            // Show modal selection editor dialog.
            // Differs from compoundEditor only in tableModel used
            SelectionEditor editor = new SelectionEditor
            (
                App.getRootFrame(),
                dictEntry_
            );
            editor.setTitle(displayName_);

            editor.show();
            // Update the cached display string.
            displayValue_ = getSelectionDisplayString();
            // Signal that editing has stopped.
            editor_.stopCellEditing();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }
    //--------------------------------------------------------------------------

    /* render selection by returning name of selected entry */
    protected String getSelectionDisplayString()
    {
        String renderedValue = "";

        try
        {
            int i = dictEntry_.selection();

            IDictionaryEntry subEntry = dictEntry_.subElements()[i];
            renderedValue = subEntry.typeDescriptor().name();
            subEntry = null;
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
        return renderedValue;
    }

    //--------------------------------------------------------------------------

    /* render selection by adding value to vector */
    protected void getSelectionValueString(Vector argListVector)
    {
        try
        {
            int i = dictEntry_.selection();

            DictionaryEntryCache dictEntryCache =
                DictionaryEntryCache.New(dictEntry_.subElements()[i]);

            dictEntryCache.toStringRaw(argListVector);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }


    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}





