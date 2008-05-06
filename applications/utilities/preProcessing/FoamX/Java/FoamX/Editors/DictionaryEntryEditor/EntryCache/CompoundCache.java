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
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntryCompoundEditor;
import FoamX.Editors.DictionaryEntryEditor.*;

public class CompoundCache
    extends DictionaryEntryCache
//    implements DictionaryEntry
{
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    /** CompoundCache constructor. */
    public CompoundCache
    (
        IDictionaryEntry dictEntry,
        TypeDescriptorCache typeDescriptor
    )
    {
        super(dictEntry, typeDescriptor);

        if (dictEntry_.subElements().length == 0)
        {
            editable_ = false;
        }
    }

    //--------------------------------------------------------------------------
    /** CompoundCache constructor. */
    public CompoundCache(IDictionaryEntry dictEntry)
    {
        // Invoke the other constructor.
        this
        (
            dictEntry,
            new TypeDescriptorCache(dictEntry.typeDescriptor(), false)
        );
    }

    //--------------------------------------------------------------------------
    /** CompoundCache constructor. */
    public CompoundCache(IDictionaryEntry dictEntry, String displayName)
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
            displayValue_ = getCompoundDisplayString(3);
        }
        return displayValue_;
    }

    //--------------------------------------------------------------------------
    /**
     * Variant of toString which appends compound entries.
     * Used to convert argument dictionary to argument string.
     */
    public void toStringRaw(Vector argListVector)
    {
        getCompoundValueString(argListVector);
    }

    //--------------------------------------------------------------------------

    // Update cached value from IDictionaryEntry
    protected void getEntryValue()
    {}

    //--------------------------------------------------------------------------

    // Update IDictionaryEntry value from cached value
    protected void setEntryValue()
    {}

    //--------------------------------------------------------------------------

    protected void initialiseEditor()
    {
        // Use compound editor.
        DictionaryEntryCompoundEditor compEdit =
            new DictionaryEntryCompoundEditor();

        // Subscribe to the edit button's action event so
        // that we can invoke the appropriate edit action for
        // the compound types.
        compEdit.addActionListener
        (
            new java.awt.event.ActionListener()
            {
                public void actionPerformed
                (
                    java.awt.event.ActionEvent evt
                )
                {
                    editButtonActionPerformed(evt);
                }
            }
        );
        editor_ = compEdit;
    }

    //--------------------------------------------------------------------------

    protected void editButtonActionPerformed(java.awt.event.ActionEvent evt)
    {
        // The edit button has been pressed. The user wants to edit this entry.
        try
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
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    /* render compound by adding first maxEntries values, separated by ',' */
    protected String getCompoundDisplayString(int maxEntries)
    {
        String vecString = "";

        try
        {
            // Loop over sub-entries and generate a string representation of
            // the vector.
            IDictionaryEntry[] subEntries = dictEntry_.subElements();
            for (int i=0; i <subEntries.length; i++)
            {
                DictionaryEntryCache dictEntryCache =
                    DictionaryEntryCache.New(subEntries[i]);

                vecString += dictEntryCache.toString();

                if (i < subEntries.length - 1 && i < maxEntries)
                {
                    // Add comma if more values to follow.
                    vecString += ", ";
                }
                else if (i == maxEntries)
                {
                    // Reached maximum number of values.
                    vecString += " ...";
                    break;
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return vecString;
    }

    //--------------------------------------------------------------------------

    /* render compound by adding all values to vector */
    protected void getCompoundValueString(Vector argListVector)
    {
        try
        {
            // Loop over sub-entries and generate a string representation of
            // the vector.
            IDictionaryEntry[] subEntries = dictEntry_.subElements();

            for (int i = 0; i < subEntries.length; i++)
            {
                DictionaryEntryCache dictEntryCache =
                    DictionaryEntryCache.New(subEntries[i]);

                dictEntryCache.toStringRaw(argListVector);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }


    //--------------------------------------------------------------------------
    /**
     *  Called by the DictionaryCellEditor's getCellEditorValue
     *  method to update the value of this item after editing.
     *  Will not be called by editors for compound types since e.g.
     *  dictionary elements cannot be added/removed in FoamX.
     */
    public boolean updateValue(Object value)
    {
        //try
        //{
        //    throw new Exception("Invalid call to updateValue.");
        //}
        //catch (Exception ex)
        //{
        //    App.handleAllExceptions(ex);
        //}
        return false;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}





