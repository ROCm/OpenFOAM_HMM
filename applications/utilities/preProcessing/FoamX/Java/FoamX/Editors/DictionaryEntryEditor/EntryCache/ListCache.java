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

public class ListCache
    extends CompoundCache
//    implements DictionaryEntry
{
    //--------------------------------------------------------------------------

    // Length of list before confirming edit.
    public static final int MAX_EDIT_LENGTH = 100;

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    /** ListCache constructor. */
    public ListCache
    (
        IDictionaryEntry dictEntry,
        TypeDescriptorCache typeDescriptor
    )
    {
        super(dictEntry, typeDescriptor);

        // Reset editable. CompoundCache set it to false for zero element lists.
        editable_ = typeDescriptor_.isEditable();
    }

    //--------------------------------------------------------------------------
    /** ListCache constructor. */
    public ListCache(IDictionaryEntry dictEntry)
    {
        // Invoke the other constructor.
        this
        (
            dictEntry,
            new TypeDescriptorCache(dictEntry.typeDescriptor(), false)
        );
    }

    //--------------------------------------------------------------------------
    /** ListCache constructor. */
    public ListCache(IDictionaryEntry dictEntry, String displayName)
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
            displayValue_ = "List...";
        }

        return displayValue_;
    }

    //--------------------------------------------------------------------------
    /**
     * Variant of toString which appends compound entries.
     * Used to convert command line arguments dictionary to argument string.
     * Note that now a List will become one argument string.
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
            // Check n elems and ask before editing.
            boolean doEdit = true;

            int nSubElements = dictEntry_.nSubElements();
            if (nSubElements > MAX_EDIT_LENGTH)
            {
                if
                (
                    JOptionPane.showConfirmDialog
                    (
                        App.getRootFrame(),
                        "List is large (" + nSubElements
                        + " elements).\n"
                        + "Are you sure you want to edit it?",
                        "Confirm Edit List",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.QUESTION_MESSAGE
                    )
                    !=
                    JOptionPane.OK_OPTION
                )
                {
                    doEdit = false;
                }
            }

            // Start editing
            if (doEdit)
            {
                // Show modal list editor dialog.
                ListEditor editor = new ListEditor
                (
                    App.getRootFrame(),
                    dictEntry_
                );
                editor.setTitle(displayName_);

                editor.show();
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

    public DictionaryEntryCache addListElement()
    {
        DictionaryEntryCache listElement = null;

        try
        {
            // Add a new element to the list.
            IDictionaryEntryHolder listElementHolder =
                new IDictionaryEntryHolder();
            dictEntry_.addElement(listElementHolder);

            // If successful, create a new DictionaryEntryCache wrapper
            // object and return it.
            listElement = DictionaryEntryCache.New(listElementHolder.value);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return listElement;
    }

    //--------------------------------------------------------------------------

    public boolean removeElement(DictionaryEntryCache listElement)
    {
        boolean bRet = false;

        try
        {
            // Remove the specified sub-type.
            dictEntry_.removeElement(listElement.getDictEntry());

            // Return success.
            bRet = true;
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
            bRet = false;
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
            bRet = false;
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





