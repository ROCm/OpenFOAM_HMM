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
package FoamX.Options;

import java.util.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;

import FoamX.App;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;

public class OptionTableModel
    extends AbstractTableModel
{
    //--------------------------------------------------------------------------

    public static final int NAME_COLUMN_INDEX  = 0;
    public static final int VALUE_COLUMN_INDEX = 1;

    protected OptionsManager optionsMgr_;        // Reference to the global options manager object.
    protected java.util.List paramList_;         // List of OptionParameter objects.
    protected java.util.Hashtable paramMap_;          // Map of OptionParameter objects.
    protected String[]             columnNames_ = { "Option", "Value" };

    //--------------------------------------------------------------------------
    /**
     * DictionaryEntryTableModel constructor for compound types.
     * @param compoundEntry Reference to the compound entry object.
     */
    public OptionTableModel(OptionsManager optionsMgr)
    {
        try
        {
            optionsMgr_ = optionsMgr;
            paramList_  = new Vector(10);
            paramMap_   = new Hashtable();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void addParameter(String paramName)
    {
        OptionParameter param = new OptionParameter(paramName);
        paramList_.add(param);
        paramMap_.put(paramName, param);
    }

    //--------------------------------------------------------------------------
    //--- Abstract AbststractTableModel Methods
    //--------------------------------------------------------------------------

    public int getRowCount()
    {
        return paramList_.size();
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
            if (rowIndex>= paramList_.size())
            {
                throw new Exception("Invalid row index");
            }

            // Get OptionParameter object for this row.
            OptionParameter param = (OptionParameter)paramList_.get(rowIndex);

            if (columnIndex == NAME_COLUMN_INDEX)
            {
                // Return the OptionParameter display name.
                obj = param.getName();
            }
            else
            {
                // Return the OptionParameter object for this entry.
                obj = param;
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
        // Value column is of type DictionaryEntry.
        if (columnIndex == VALUE_COLUMN_INDEX)
        {
            return OptionParameter.class;
        }

        // Name column is of type String.
        return String.class;
    }

    //--------------------------------------------------------------------------

    public boolean isCellEditable(int rowIndex, int columnIndex)
    {
        boolean bEditable = false;

        try
        {
            if (columnIndex == VALUE_COLUMN_INDEX && rowIndex <paramList_.size())
            {
                // Get OptionParameter object for this row.
                OptionParameter param = (OptionParameter)paramList_.get(rowIndex);

                // Return the editable flag.
                bEditable = param.isEditable();
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
        // Ignore any incoming DictionaryEntry objects.
        fireTableCellUpdated(rowIndex, columnIndex);
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}


