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
import javax.swing.tree.*;

import FoamX.App;

public class FoamXTreeRenderer
    extends JLabel
    implements TreeCellRenderer
{
    //-----------------------------------------------------------------------------------------

    private Icon defaultIcon_;

    //-----------------------------------------------------------------------------------------

    /** FoamXTreeRenderer constructor. */
    public FoamXTreeRenderer()
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

    public Component getTreeCellRendererComponent
    (
        JTree   tree,
        Object  value,
        boolean selected,
        boolean expanded,
        boolean leaf,
        int     row,
        boolean hasFocus
    )
    {
        try
        {
            // User object should support the FoamXTreeItem interface.
            DefaultMutableTreeNode nodeItem = (DefaultMutableTreeNode)value;

            // See if this object has a user object associated with it.
            if (nodeItem.getUserObject() != null)
            {
                FoamXTreeItem item = (FoamXTreeItem)nodeItem.getUserObject();

                // Set label text and icon.
                setText(item.getText());
                setIcon(item.getIcon());

                if (selected)
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
                    setBackground(tree.getBackground());
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return this;
    }

    //-----------------------------------------------------------------------------------------

    public void paint(Graphics g)
    {
        int textstart = 0;

        if (getIcon() != null)
        {
            textstart = getIcon().getIconWidth() + getIconTextGap();
            getIcon().paintIcon(this, g, 0, 0);
        }

        int textwidth = getSize().width - textstart;
        g.setColor(getBackground());
        g.fillRect(textstart, 0, textwidth, getSize().height);
        g.setColor(getForeground());
        g.drawString(getText(), textstart, (getSize().height + g.getFontMetrics().getAscent()) / 2);
    }

    //-----------------------------------------------------------------------------------------
    /** This interface defines the methods required by all items using the FoamXTreeRenderer
     *  tree renderer class.
     */
    public static interface FoamXTreeItem
    {
        public Icon getIcon();
        public String getText();
    }

    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
}





