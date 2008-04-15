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
import FoamX.Util.FoamXAny;
import FoamX.Util.FoamXTypeEx;
import FoamX.Editors.TypeEditor.TypeDescriptorCache;

public class StringCache
    extends DictionaryEntryCache
//    implements DictionaryEntry
{

    //--------------------------------------------------------------------------
    /** StringCache constructor. */
    public StringCache
    (
        IDictionaryEntry dictEntry,
        TypeDescriptorCache typeDescriptor
    )
    {
        super(dictEntry, typeDescriptor);
    }

    //--------------------------------------------------------------------------
    /** StringCache constructor. */
    public StringCache(IDictionaryEntry dictEntry)
    {
        // Invoke the other constructor.
        this
        (
            dictEntry,
            new TypeDescriptorCache(dictEntry.typeDescriptor(), false)
        );
    }

    //--------------------------------------------------------------------------
    /** StringCache constructor. */
    public StringCache(IDictionaryEntry dictEntry, String displayName)
    {
        // Invoke the other constructor.
        this(dictEntry);

        // Use the given display name.
        displayName_ = displayName;
    }

    //--------------------------------------------------------------------------

    protected void initialiseEditor()
    {
        // Give dictionaryEntryCache opportunity to install combo editor
        super.initialiseEditor();

        if (editor_ == null)
        {
            // No custom editor - use string-based
            // DefaultCellEditor.
        }
    }

    //--------------------------------------------------------------------------
    //---- DictionaryEntry Interface Methods
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    /**
     *  Called by the DictionaryCellEditor's getCellEditorValue
     *  method to update the value of this item after editing.
     *  Will not be called by editors for compound types since e.g.
     *  dictionary elements cannot be added/removed in FoamX.
     */
    public boolean updateValue(Object value)
    {
        boolean bRet = false;
        String strValue;

        try
        {
            if (value != null)
            {
                // Incoming value object is a string.
                strValue = (String)value;

                // Update the Any object.
                value_.setValue(strValue);
                // Send the Any object to the dictionary entry object.
                setEntryValue();

                bRet = true;
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return bRet;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}





