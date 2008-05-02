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
import FoamX.Editors.FixedListEditor;
import FoamX.Editors.TypeEditor.TypeDescriptorCache;

public class FixedListCache
    extends CompoundCache
//    implements DictionaryEntry
{
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    /** FixedListCache constructor. */
    public FixedListCache
    (
        IDictionaryEntry dictEntry,
        TypeDescriptorCache typeDescriptor
    )
    {
        super(dictEntry, typeDescriptor);
    }

    //--------------------------------------------------------------------------
    /** FixedListCache constructor. */
    public FixedListCache(IDictionaryEntry dictEntry)
    {
        // Invoke the other constructor.
        this
        (
            dictEntry,
            new TypeDescriptorCache(dictEntry.typeDescriptor(), false)
        );
    }

    //--------------------------------------------------------------------------
    /** FixedListCache constructor. */
    public FixedListCache(IDictionaryEntry dictEntry, String displayName)
    {
        // Invoke the other constructor.
        this(dictEntry);

        // Use the given display name.
        displayName_ = displayName;
    }

    //--------------------------------------------------------------------------
    /**
     * Variant of toString which appends compound entries.
     * Used to convert command line arguments dictionary to argument string.
     * Note that now a FixedList will become one argument string.
     */
    public void toStringRaw(Vector argListVector)
    {
        // Collect arguments into vector (only so we can reuse
        // CompoundCache.toStringRaw functionality)
        Vector subArgs = new Vector();
        super.toStringRaw(subArgs);

        // Compose single string out of arguments. Enclose in brackets.

        String argString = "(";

        for (int i = 0; i < subArgs.size(); i++)
        {
            if (i > 0)
            {
                argString += " ";
            }
            argString += (String)subArgs.elementAt(i);
        }

        argString += ")";

        // Add as single argument
        argListVector.addElement(argString);
    }
    
    //--------------------------------------------------------------------------

    protected void editButtonActionPerformed(java.awt.event.ActionEvent evt)
    {
        // The edit button has been pressed. The user wants to edit this entry.
        try
        {
            // Show modal vector space editor dialog.
            FixedListEditor editor = new FixedListEditor
            (
                App.getRootFrame(),
                dictEntry_,
                typeDescriptor_
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
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}





