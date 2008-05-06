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

import java.awt.event.ActionEvent;
import javax.swing.*;
import javax.swing.tree.*;
import javax.swing.event.*;
import java.util.Enumeration;

import FoamX.App;

public class ExpandTreeAction
    extends AbstractAction
{
    //--------------------------------------------------------------------------

    protected JTree tree_;

    //--------------------------------------------------------------------------
    /** ExpandTreeAction constructor. */
    public ExpandTreeAction(JTree tree)
    {
        tree_ = tree;

        putValue(Action.SMALL_ICON, App.getResources().getIcon("Tree.ExpandTreeImage"));
        putValue(Action.NAME, "ExpandAll");
        putValue(Action.SHORT_DESCRIPTION, "Expand All");
        putValue(Action.LONG_DESCRIPTION, "Expand All");
    }

    //--------------------------------------------------------------------------

    public void actionPerformed(ActionEvent e)
    {
        try
        {
            // Get tree model.
            TreeModel model = tree_.getModel();

            DefaultMutableTreeNode root = (DefaultMutableTreeNode)model.getRoot();
            Enumeration iter = root.depthFirstEnumeration();
            while (iter.hasMoreElements())
            {
                DefaultMutableTreeNode node = (DefaultMutableTreeNode)iter.nextElement();
                TreePath tp = new TreePath(node.getPath());
                tree_.expandPath(tp);
            }
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
