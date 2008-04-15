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
package FoamX.Editors.TypeEditor;

import java.awt.Component;
import javax.swing.*;
import javax.swing.table.*;

class TypeParameterEditor
    extends DefaultCellEditor
{
    //--------------------------------------------------------------------------

    private TypeParameter   typeParam_    = null;
    private TableCellEditor customEditor_ = null;

    //--------------------------------------------------------------------------
    /** TypeParameterEditor constructor. */
    public TypeParameterEditor()
    {
        super(new JTextField());
    }

    //--------------------------------------------------------------------------
    /** Implementation of TableCellEditor interface's getTableCellEditorComponent function. */

    public Component getTableCellEditorComponent
    (
        JTable  table,
        Object  value,
        boolean isSelected,
        int     row,
        int     column
    )
    {
        typeParam_ = (TypeParameter)value;

        // Delegate to custom editor if required.
        if (typeParam_.customEditor() != null)
        {
            customEditor_ = typeParam_.customEditor();
            return customEditor_.getTableCellEditorComponent(table, value, isSelected, row, column);
        }
        else if (typeParam_.entryEditor() != null)
        {
            customEditor_ = typeParam_.entryEditor();
            return customEditor_.getTableCellEditorComponent
            (
                table, typeParam_.getValue(), isSelected, row, column
            );
        }

        // Otherwise, return the default editor component and edit the value string.
        return super.getTableCellEditorComponent
        (
            table, typeParam_.getValueString(), isSelected, row, column
        );
    }

    //--------------------------------------------------------------------------

    public Object getCellEditorValue()
    {
        if (customEditor_ != null)
        {
            // Delegate to custom editor.
            Object objValue = customEditor_.getCellEditorValue();
        }
        else
        {
            // Intercept the "value" object from this editor.
            // Of class String for DefaultCellEditor (JTextField) class.
            String objValue = (String)super.getCellEditorValue();

            // Update the TypeParameter object's value.
            typeParam_.updateValue(objValue);
        }

        // Return the TypeParameter object being edited.
        // This is the value that gets sent to the model.
        return typeParam_;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}


