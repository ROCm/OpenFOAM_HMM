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

package FoamX.Editors.ApplicationEditor;

import javax.swing.*;
import javax.swing.tree.*;
import java.awt.Component;
import FoamX.App;

public class BoundaryDefinitionTreeRenderer
    extends FoamX.Util.FoamXTreeRenderer
{
    // Icons.
    private Icon[] icons_;

    // Symbolic constants for icon indices.
    public static final int UNDEFINED_ICON   = 0;
    public static final int ROOT_ICON        = 1;
    public static final int PATCHTYPE_ICON   = 2;
    public static final int BOUNDARYDEF_ICON = 3;

    //-----------------------------------------------------------------------------------------

    public BoundaryDefinitionTreeRenderer()
    {
        try
        {
            // Load icons.
            icons_ = new Icon[4];
            icons_[UNDEFINED_ICON]   = App.getResources().getIcon("UndefinedImage");
            icons_[ROOT_ICON]        = App.getResources().getIcon("PatchRootImage");
            icons_[PATCHTYPE_ICON]   = App.getResources().getIcon("PatchTypeImage");
            icons_[BOUNDARYDEF_ICON] = App.getResources().getIcon("BoundaryDefinitionImage");
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
            super.getTreeCellRendererComponent(tree,
                                                value,
                                                selected,
                                                expanded,
                                                leaf,
                                                row,
                                                hasFocus);

            DefaultMutableTreeNode      treeNode = (DefaultMutableTreeNode)value;
            BoundaryDefinitionModelItem nodeInfo = (BoundaryDefinitionModelItem)treeNode.getUserObject();

            // Set label text.
            setText(nodeInfo.toString());

            // Set an appropriate icon according to the type of node.
            if (nodeInfo.isRoot())
            {
                setIcon(icons_[ROOT_ICON]);
            }
            else if (nodeInfo.isPatchType())
            {
                setIcon(icons_[PATCHTYPE_ICON]);
            }
            else if (nodeInfo.isBoundaryDefinition())
            {
                setIcon(icons_[BOUNDARYDEF_ICON]);
            }
            else
            {
                setIcon(icons_[UNDEFINED_ICON]);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return this;
    }

    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------
}




