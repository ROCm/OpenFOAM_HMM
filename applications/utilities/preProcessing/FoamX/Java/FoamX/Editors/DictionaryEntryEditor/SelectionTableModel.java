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
package FoamX.Editors.DictionaryEntryEditor;

import java.util.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;

import FoamX.App;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;
import FoamXServer.IDictionaryEntry;
import FoamX.Editors.DictionaryEntryEditor.EntryCache.DictionaryEntryCache;
public class SelectionTableModel
    extends AbstractTableModel
    implements CompoundEntryListener
{
    //--------------------------------------------------------------------------

    public static final int NAME_COLUMN_INDEX  = 0;
    public static final int VALUE_COLUMN_INDEX = 1;

    // Root DictionaryEntry object for this table.
    protected DictionaryEntryCache rootEntry_;
    protected String dictionaryKey_;

    // List of DictionaryEntry objects.
    protected java.util.List dictEntryList_;

    // Map of DictionaryEntry objects.
    protected java.util.Hashtable dictEntryMap_;
    protected EventListenerList listenerList_;
    protected String[] columnNames_ = { "Name", "Value"};

    //--------------------------------------------------------------------------
    /** SelectionTableModel default constructor for derived classes..
     */
    protected SelectionTableModel()
    {
        try
        {
            rootEntry_     = null;
            dictionaryKey_ = "";
            dictEntryList_ = new Vector(10);
            dictEntryMap_  = new Hashtable();
            listenerList_  = new EventListenerList();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /**
     * SelectionTableModel constructor for compound types.
     * @param compoundEntry Reference to the compound entry object.
     */
    public SelectionTableModel(IDictionaryEntry compoundEntry)
    {
        // Default construction.
        this();

        try
        {
            // Store reference to IDictionaryEntry for this compound entry.
            rootEntry_     = DictionaryEntryCache.New(compoundEntry);
            // Default dictionary key.
            dictionaryKey_ = rootEntry_.getTypeDescriptor().getName();

            // Initialise event listener list.
            listenerList_ = new EventListenerList();

            // Create a DictionaryEntryCache object for each sub-entry.
            IDictionaryEntry[] subEntries = compoundEntry.subElements();

            dictEntryList_ = new Vector(subEntries.length);
            for (int i=0; i <subEntries.length; i++)
            {
                // Create DictionaryEntryCache object.
                IDictionaryEntry subEntry = subEntries[i];
                DictionaryEntry dictEntry = DictionaryEntryCache.New(subEntry);

                // Add to list.
                dictEntryList_.add(i, dictEntry);

                // Add to map.
                dictEntryMap_.put(dictEntry.getEntryName(), dictEntry);

                // Subscribe to the openSubDictionary event so that we can
                // forward the event to the listeners of this object.
                dictEntry.addCompoundEntryListener(this);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public DictionaryEntryCache getRootEntry()
    {
        return rootEntry_;
    }

    public String getNameColumnTitle()
    {
        return columnNames_[NAME_COLUMN_INDEX];
    }
    public void setNameColumnTitle(String title)
    {
        columnNames_[NAME_COLUMN_INDEX] = title;
    }

    public String getValueColumnTitle()
    {
        return columnNames_[VALUE_COLUMN_INDEX];
    }
    public void setValueColumnTitle(String title)
    {
        columnNames_[VALUE_COLUMN_INDEX] = title;
    }

    public String getDictionaryKey()
    {
        return dictionaryKey_;
    }
    public void setDictionaryKey(String key)
    {
        dictionaryKey_ = key;
    }

    //--------------------------------------------------------------------------
    /** Returns the index of the specified entry. */
    public int getIndexOfEntry(String entryName)
    {
        // Return zero by default.
        int index = 0;

        for (int i = 0; i <dictEntryList_.size(); i++)
        {
            DictionaryEntryCache dictEntry =
                (DictionaryEntryCache)dictEntryList_.get(i);

            if (entryName.compareToIgnoreCase(dictEntry.getEntryName()) == 0)
            {
                index = i;
                break;
            }
        }

        return index;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
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

    protected boolean fireOpenSubDictionaryEvent(CompoundEntryEvent evt)
    {
        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();

        // Prepend the dictionary key and send to the module.
        String fullDictName = dictionaryKey_ + "/" + evt.getName();

        // Create the event object.
        CompoundEntryEvent event =
            new CompoundEntryEvent
            (
                this,
                CompoundEntryEvent.TYPE_DICTIONARY,
                fullDictName
            );

        // Fire the event.
        if (evt != null)
        {
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
        }
        else
        {
            App.printMessage
            (
                "SelectionTableModel.fireOpenSubDictionaryEvent"
                + "(CompoundEntryEvent evt)"
                + " : cannot fire event"
            );
        }

        return false;
    }

    //--------------------------------------------------------------------------
    //--- CompoundEntryListener Interface Methods
    //--------------------------------------------------------------------------

    public boolean openSubDictionary(CompoundEntryEvent evt)
    {
        // Bounce the openSubDictionary event to any listeners.
        return fireOpenSubDictionaryEvent(evt);
    }

    //--------------------------------------------------------------------------
    //--- Abstract AbststractTableModel Methods
    //--------------------------------------------------------------------------

    public int getRowCount()
    {
        return dictEntryList_.size();
    }

    //--------------------------------------------------------------------------

    public int getColumnCount()
    {
        return columnNames_.length;
    }

    //--------------------------------------------------------------------------

    public Object getValueAt(int rowIndex, int columnIndex)
    {
        Object obj = null;

        try
        {
            // Check for valid row index.
            if (rowIndex>= dictEntryList_.size())
            {
                throw new Exception("Invalid row index");
            }

            // Get DictionaryEntry object for this row.
            DictionaryEntryCache dictEntry =
                (DictionaryEntryCache)dictEntryList_.get(rowIndex);

            if (columnIndex == NAME_COLUMN_INDEX)
            {
                // Return the DictionaryEntry display name.
                obj = dictEntry.getEntryName();
            }
            else if (columnIndex == VALUE_COLUMN_INDEX)
            {
                obj = dictEntry;
            }
            else
            {
                throw new Exception("Invalid column index");
            }

        }
        catch (FoamXIOError ioErr)
        {
            App.handleException(ioErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return obj;
    }

    //--------------------------------------------------------------------------
    //--- Overriden AbststractTableModel Methods
    //--------------------------------------------------------------------------

    public String getColumnName(int columnIndex)
    {
        return columnNames_[columnIndex];
    }

    //--------------------------------------------------------------------------

    public Class getColumnClass(int columnIndex)
    {
        // Render according to value of row 0 ( so 0=string, 1=bool)
        return getValueAt(0, columnIndex).getClass();
    }

    //--------------------------------------------------------------------------

    public boolean isCellEditable(int rowIndex, int columnIndex)
    {
        boolean bEditable = false;

        try
        {
            if (columnIndex == VALUE_COLUMN_INDEX)
            {
                // Get DictionaryEntry object for this row.
                DictionaryEntry dictEntry =
                    (DictionaryEntry)dictEntryList_.get(rowIndex);

                // Return the editable flag.
                bEditable = dictEntry.isEditable();
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return bEditable;
    }

    //--------------------------------------------------------------------------

    public void setValueAt(Object aValue, int rowIndex, int columnIndex)
    {
        try
        {
            if (columnIndex == NAME_COLUMN_INDEX)
            {
            }
            else if (columnIndex == VALUE_COLUMN_INDEX)
            {
                // No need to change current selection here
                // (Current-selection set by selection changes, not by editing)

                // Ignore any incoming DictionaryEntry objects.
                fireTableCellUpdated(rowIndex, columnIndex);
            }
            else
            {
                throw new Exception("Invalid column index");
            }
        }
        catch (FoamXIOError ioErr)
        {
            App.handleException(ioErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
