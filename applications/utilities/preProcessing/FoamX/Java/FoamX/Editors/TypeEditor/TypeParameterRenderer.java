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

import javax.swing.JTable;
import javax.swing.table.*;
import java.awt.Component;

class TypeParameterRenderer
    extends DefaultTableCellRenderer
{

    //--------------------------------------------------------------------------

    /** TypeParameterRenderer constructor. */
    public TypeParameterRenderer()
    {}

    //--------------------------------------------------------------------------

    public Component getTableCellRendererComponent
    (
        JTable  table,
        Object  value,
        boolean isSelected,
        boolean hasFocus,
        int     row,
        int     column
    )
    {
        TypeParameter typeParam = (TypeParameter)value;

        // Use the value to render this cell.
        String    valueText = typeParam.getValueString();
        Component comp      = super.getTableCellRendererComponent(table, valueText, isSelected, hasFocus, row, column);

        // If this parameter is read-only, grey it out.
        if (!typeParam.isEditable())
        {
            comp.setForeground(java.awt.Color.gray);
        } else {
            comp.setForeground(java.awt.Color.black);
        }

        return comp;
    }

    //--------------------------------------------------------------------------
}



