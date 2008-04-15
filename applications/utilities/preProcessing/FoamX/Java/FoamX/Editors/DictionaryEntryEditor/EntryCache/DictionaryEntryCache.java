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
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntry;
import FoamX.Editors.DictionaryEntryEditor.CompoundEntryListener;
import FoamX.Editors.DictionaryEntryEditor.CompoundEntryEvent;
import FoamX.Editors.DictionaryEntryEditor.DictionaryEntryCompoundEditor;

public class DictionaryEntryCache
    implements DictionaryEntry
{
    //--------------------------------------------------------------------------

    // Reference to the dict entry object for this item.
    protected IDictionaryEntry dictEntry_;

    // Cached type descriptor for this entry.
    protected TypeDescriptorCache typeDescriptor_;
    // Cached value
    protected FoamXAny value_;
    // Cached type of type descriptor
    protected int typeValue_;

    // Name of this entry. Defaults to the display name specified by the
    // type descriptor.
    protected String displayName_;

    // Cached String representation of the current value of this entry.
    protected String displayValue_;

    protected boolean editable_;
    protected boolean current_;
    protected boolean modified_;

    // Custom renderer object.
    protected TableCellRenderer renderer_;

    // Custom editor object.
    protected TableCellEditor editor_;

    // Listeners for open sub dictionary events
    protected EventListenerList listenerList_;


    //--------------------------------------------------------------------------

    /** DictionaryEntryCache constructor. */
    public DictionaryEntryCache
    (
        IDictionaryEntry dictEntry,
        TypeDescriptorCache typeDescriptor
    )
    {
        try
        {
            // Store IDictionaryEntry reference.
            dictEntry_ = dictEntry;

            // Get TypeDescriptor information.
            typeDescriptor_ = typeDescriptor;
            FoamXType type  = typeDescriptor_.getType();
            typeValue_      = type.value();
            displayName_    = typeDescriptor_.getDisplayName();
            displayValue_   = "";

            // Make sure that we have a name to display.
            if (displayName_ == null || displayName_.length() == 0)
            {
                displayName_ = typeDescriptor_.getName();
            }

            // Initialise event listener list.
            listenerList_ = new EventListenerList();

            // Get current value if we have one. Need a value before we
            // initialise the Renderer and Editor objects.
            getEntryValue();

            // Editable if typeDescriptor says so
            editable_ = typeDescriptor_.isEditable();

            // Now, initialise rest of member variables.
            current_ = false;
            modified_ = dictEntry.modified();

            // Initialise the renderer and editor for this entry.
            initialiseRenderer();
            initialiseEditor();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** DictionaryEntryCache constructor. */
    public DictionaryEntryCache(IDictionaryEntry dictEntry)
    {
        // Invoke the other constructor.
        this
        (
            dictEntry,
            new TypeDescriptorCache(dictEntry.typeDescriptor(), false)
        );
    }

    //--------------------------------------------------------------------------
    /** DictionaryEntryCache constructor. */
    public DictionaryEntryCache(IDictionaryEntry dictEntry, String displayName)
    {
        // Invoke the other constructor.
        this(dictEntry);

        // Use the given display name.
        displayName_ = displayName;
    }

    //--------------------------------------------------------------------------

    /** DictionaryEntryCache factory. */
    static public DictionaryEntryCache New
    (
        IDictionaryEntry dictEntry,
        TypeDescriptorCache typeDescriptor
    )
    {
        try
        {
            int typeValue = typeDescriptor.getType().value();

            switch (typeValue)
            {
                case FoamXType._Type_Boolean:
                    return new BooleanCache(dictEntry, typeDescriptor);
                case FoamXType._Type_Label:
                    return new LabelCache(dictEntry, typeDescriptor);
                case FoamXType._Type_Scalar:
                    return new ScalarCache(dictEntry, typeDescriptor);
                case FoamXType._Type_Char:
                    return new CharCache(dictEntry, typeDescriptor);
                case FoamXType._Type_Word:
                    return new WordCache(dictEntry, typeDescriptor);
                case FoamXType._Type_String:
                    return new StringCache(dictEntry, typeDescriptor);

                case FoamXType._Type_RootDir:
                    return new RootDirCache(dictEntry, typeDescriptor);
                case FoamXType._Type_RootAndCase:
                    return new RootAndCaseCache(dictEntry, typeDescriptor);
                case FoamXType._Type_CaseName:
                    return new CaseNameCache(dictEntry, typeDescriptor);
                case FoamXType._Type_HostName:
                    return new HostNameCache(dictEntry, typeDescriptor);
                case FoamXType._Type_File:
                    return new StringCache(dictEntry, typeDescriptor);
                case FoamXType._Type_Directory:
                    return new StringCache(dictEntry, typeDescriptor);
                case FoamXType._Type_Time:
                    return new TimeCache(dictEntry, typeDescriptor);

                case FoamXType._Type_DimensionSet:
                    return new DimensionSetCache(dictEntry, typeDescriptor);
                case FoamXType._Type_FixedList:
                    return new FixedListCache(dictEntry, typeDescriptor);
                case FoamXType._Type_List:
                    return new ListCache(dictEntry, typeDescriptor);
                case FoamXType._Type_Selection:
                    if (dictEntry.subElements().length == 1)
                    {
                        // Selection with one entry.
                        App.printMessage
                        (
                            App.DEBUGLEVEL_DEBUG,
                            "Skipping single element Selection " + dictEntry
                        );

                        return New
                        (
                            dictEntry.subElements()[0],
                            new TypeDescriptorCache
                            (
                                typeDescriptor.TypeDescriptor().subTypes()[0],
                                false
                            )
                        );
                    }
                    else
                    {
                        return new SelectionCache(dictEntry, typeDescriptor);
                    }
                case FoamXType._Type_Dictionary:
                    return new DictionaryCache(dictEntry, typeDescriptor);
                case FoamXType._Type_Compound:
                    return new CompoundCache(dictEntry, typeDescriptor);
                default:
                    throw new Exception
                    (
                        "Invalid type in DictionaryEntryCache::New."
                    );
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
        return null;
    }

    //--------------------------------------------------------------------------
    /** DictionaryEntryCache factory. */
    static public DictionaryEntryCache New(IDictionaryEntry dictEntry)
    {
        // Invoke the other factory.
        return New
        (
            dictEntry,
            new TypeDescriptorCache(dictEntry.typeDescriptor(), false)
        );
    }

    //--------------------------------------------------------------------------
    /** DictionaryEntryCache constructor. */
    static public DictionaryEntryCache New
    (
        IDictionaryEntry dictEntry,
        String displayName
    )
    {
        // Invoke the other constructor.
        DictionaryEntryCache entry = New(dictEntry);

        // Use the given display name.
        entry.setEntryName(displayName);

        return entry;
    }

    //--------------------------------------------------------------------------

    public IDictionaryEntry getDictEntry()
    {
        return dictEntry_;
    }
    public TypeDescriptorCache getTypeDescriptor()
    {
        return typeDescriptor_;
    }
    public void setEditable(boolean editable)
    {
        editable_ = editable;
    }
    public void setCurrent(boolean current)
    {
        current_ = current;
    }
    public boolean isModified()
    {
        return modified_;
    }

    public void setEntryName(String displayName)
    {
        displayName_ = displayName;
    }

    //--------------------------------------------------------------------------
    /**
     * Get array of sub elements
     */
    public DictionaryEntryCache[] getSubElements()
    {
        IDictionaryEntry[] subEntries = dictEntry_.subElements();

        DictionaryEntryCache[] cachedElems =
            new DictionaryEntryCache[subEntries.length];

        for (int i=0; i<subEntries.length; i++)
        {
            cachedElems[i] = DictionaryEntryCache.New(subEntries[i]);
        }

        return cachedElems;
    }

    //--------------------------------------------------------------------------
    /**
     * Return subdictionary by name
     */
    public DictionaryEntryCache getSubEntry
    (
        String subName
    )
    {
        IDictionaryEntry[] subEntries = dictEntry_.subElements();
        for (int i=0; i <subEntries.length; i++)
        {
            IDictionaryEntry subEntry = subEntries[i];
            if (subEntry.typeDescriptor().name().equals(subName))
            {
                // Release Corba object
                subEntries = null;

                return DictionaryEntryCache.New(subEntry);
            }
        }

        // Release Corba object
        subEntries = null;

        return null;
    }
    
    //--------------------------------------------------------------------------
    /**
     * Write for debugging.
     */
    public void write()
    {
        try
        {
            System.out.println("displayName_:" + displayName_);
            System.out.println("value_:" + value_);
            System.out.println("typeValue_:" + typeValue_);
            System.out.println("displayValue_:" + displayValue_);
            System.out.println("editable_:" + editable_);
            System.out.println("current_:" + current_);
            System.out.println("modified_:" + modified_);
            //System.out.println("dictEntry_:" + dictEntry_);

            System.out.println("typeDescriptor_:" + typeDescriptor_);
            String[] valueList = typeDescriptor_.getValueList();
            if (valueList != null)
            {
                for (int i=0; i <valueList.length; i++)
                {
                    System.out.println("    value:" + valueList[i]);
                }
            }

            IDictionaryEntry[] subEntries = dictEntry_.subElements();
            for (int i=0; i <subEntries.length; i++)
            {
                DictionaryEntryCache subEntry = new DictionaryEntryCache
                (
                    subEntries[i]
                );
                System.out.println("    subEntry:" + subEntry);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
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

    // Fill cached value
    protected void getEntryValue()
    {
        value_ = new FoamXAny(dictEntry_.value());
    }

    //--------------------------------------------------------------------------

    // Write cached value
    protected void setEntryValue()
    {
        try
        {
            dictEntry_.value(value_.getAny());

            modified_ = true;
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    protected void initialiseRenderer()
    {
        // No custom renderer - use DefaultTableCellRenderer wich
        // renders the name of this object, via the toString method,
        // into a JLabel object.
    }

    //--------------------------------------------------------------------------

    protected void initialiseEditor()
    {
        try
        {
            // If this type has a constrained value list, use a combo box
            // to edit it.
            if (typeDescriptor_.getValueList().length > 0)
            {
                String[] valueList = typeDescriptor_.getValueList();
                JComboBox combo = new JComboBox();
                combo.setFont(new java.awt.Font("Dialog", 0, 10));
                for (int i=0; i <valueList.length; i++)
                {
                    combo.addItem(valueList[i]);
                }
                editor_ = new DefaultCellEditor(combo);
            }
            else
            {
                // No custom editor - use string-based
                // DefaultCellEditor.
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }


    //--------------------------------------------------------------------------
    //---- DictionaryEntry Interface Methods
    //--------------------------------------------------------------------------

    public String getEntryName()
    {
        return displayName_;
    }

    public String getEntryDescription()
    {
        return typeDescriptor_.getDescription();
    }

    public String getEntryValueString()
    {
        return toString();
    }

    public boolean isEditable()
    {
        return editable_;
    }

    public boolean isCurrent()
    {
        return current_;
    }

    public TableCellRenderer getRenderer()
    {
        return renderer_;
    }

    public TableCellEditor getEditor()
    {
        return editor_;
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
        try
        {
             throw new Exception("updateValue not overloaded!");
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
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
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}





