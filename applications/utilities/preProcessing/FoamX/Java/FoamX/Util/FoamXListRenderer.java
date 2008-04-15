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
package FoamX.Util;

import java.awt.*;
import javax.swing.*;

import FoamX.App;

public class FoamXListRenderer
    extends DefaultListCellRenderer
{
    //-----------------------------------------------------------------------------------------

    private Icon defaultIcon_;

    //-----------------------------------------------------------------------------------------
    /** FoamXListRenderer constrcutor. */
    public FoamXListRenderer()
    {
        try
        {
            // Get a reference to the main resource helper and load a default icon.
            defaultIcon_ = new ImageIcon(App.getResources().getResourceURL("UndefinedImage"));
            setIcon(defaultIcon_);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //-----------------------------------------------------------------------------------------

    public Component getListCellRendererComponent
    (
        JList   list,
        Object  value,
        int     index,
        boolean isSelected,
        boolean cellHasFocus
    )
    {
        try
        {
            Component retValue = super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);

            // Object should support the FoamXListItem interface.
            FoamXListItem item = (FoamXListItem)value;

            // Set label text and icon.
            setText(item.getText());
            setIcon(item.getIcon());

            if (isSelected)
            {
                // Hack for X where the system colors are not all defined properly.
                //if (SystemColor.textHighlightText.getRGB() == SystemColor.textHighlight.getRGB())
                //{
                //    setForeground(Color.white);
                //    setBackground(Color.blue);
                //}
                //else
                {
                    setForeground(SystemColor.textHighlightText);
                    setBackground(SystemColor.textHighlight);
                }
            }
            else
            {
                // Hack for X where the system colors are not all defined properly.
                //if (SystemColor.textText.getRGB() == tree.getBackground().getRGB())
                //{
                //    setForeground(Color.black);
                //}
                //else
                {
                    setForeground(SystemColor.textText);
                }
                setBackground(list.getBackground());
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return this;
    }

    //-----------------------------------------------------------------------------------------
    /** This interface defines the methods required by all items using the FoamXListRenderer
     *  list renderer class.
     */
    public static interface FoamXListItem
    {
        public Icon getIcon();
        public String getText();
    }

    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
}


