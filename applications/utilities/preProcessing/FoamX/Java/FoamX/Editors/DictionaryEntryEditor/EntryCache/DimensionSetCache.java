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
import FoamX.Editors.DimensionSetEditor;
import FoamX.Util.FoamXAny;
import FoamX.Util.FoamXTypeEx;
import FoamX.Editors.TypeEditor.TypeDescriptorCache;

public class DimensionSetCache
    extends DictionaryEntryCache
{
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    /** DimensionSetCache constructor. */
    public DimensionSetCache
    (
        IDictionaryEntry dictEntry,
        TypeDescriptorCache typeDescriptor
    )
    {
        super(dictEntry, typeDescriptor);
    }

    //--------------------------------------------------------------------------
    /** DimensionSetCache constructor. */
    public DimensionSetCache(IDictionaryEntry dictEntry)
    {
        // Invoke the other constructor.
        this
        (
            dictEntry,
            new TypeDescriptorCache(dictEntry.typeDescriptor(), false)
        );
    }

    //--------------------------------------------------------------------------
    /** DimensionSetCache constructor. */
    public DimensionSetCache(IDictionaryEntry dictEntry, String displayName)
    {
        // Invoke the other constructor.
        this(dictEntry);

        // Use the given display name.
        displayName_ = displayName;
    }
    

    //--------------------------------------------------------------------------
    /**
     * Return an appropriate string representation of the
     * current value in the unselected state. Used by the
     * default cell renderer to set the text for the label.
     */
    public String toString()
    {
        try
        {
            // Convert the Any object to a nice string.
            displayValue_ = value_.toString();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
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
        try
        {
            // Convert the Any object to a nice string.
            argListVector.addElement(value_.toString());
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    protected void editButtonActionPerformed(java.awt.event.ActionEvent evt)
    {
        // The edit button has been pressed. The user wants to edit this entry.
        try
        {
            // Show modal dimension set editor dialog.
            DimensionSetEditor editor =
                new DimensionSetEditor(App.getRootFrame(), dictEntry_);

            editor.setTitle(displayName_);

            editor.show();
            // Update value in case it has been updated.
            getEntryValue();
            // Signal that editing has stopped.
            editor_.stopCellEditing();
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
    //--------------------------------------------------------------------------
}





