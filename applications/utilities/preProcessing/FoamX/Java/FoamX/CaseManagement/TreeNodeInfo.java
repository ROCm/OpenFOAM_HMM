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
package FoamX.CaseManagement;

import java.awt.event.*;
import javax.swing.*;
import javax.swing.tree.DefaultMutableTreeNode;

import FoamX.App;
import FoamX.Util.MenuItemActionChangedListener;

class TreeNodeInfo
    implements FoamX.Util.FoamXTreeRenderer.FoamXTreeItem
{
    //--------------------------------------------------------------------------

    protected String nodeLabel_;
    protected int nodeID_;
    protected JPopupMenu popupMenu_;
    protected Action defaultAction_;
    protected String defaultContext_;
    protected Icon icon_;

    //--------------------------------------------------------------------------
    /** ModuleTreeNode constructor. */
    TreeNodeInfo(String nodeLabel, int nodeID, Icon icon)
    {
        nodeLabel_     = nodeLabel;
        nodeID_        = nodeID;
        popupMenu_     = null;
        defaultAction_ = null;
        icon_          = icon;
    }

    //--------------------------------------------------------------------------

    public String toString()
    {
        return nodeLabel_;
    }

    String getNodeLabel()
    {
        return nodeLabel_;
    }
    void setNodeLabel(String label)
    {
        nodeLabel_ = label;
    }

    int getNodeID()
    {
        return nodeID_;
    }

    JPopupMenu getContextMenu()
    {
        return popupMenu_;
    }
    void setContextMenu(JPopupMenu contextMenu)
    {
        popupMenu_ = contextMenu;
    }

    //--------------------------------------------------------------------------

    void addAction(Action action, String context)
    {
        try
        {
            // Create the popup menu if required.
            if (popupMenu_ == null)
            {
                popupMenu_ = new JPopupMenu();
                popupMenu_.setFont(new java.awt.Font("Dialog", 0, 10));

                // First action is the default action.
                defaultAction_  = action;
                defaultContext_ = context;
            }

            // Add new menu item.
            JMenuItem menuItem = popupMenu_.add(action);
            menuItem.setFont(popupMenu_.getFont());

            // Set action command and link up action listener object.
            menuItem.setActionCommand(context);
            action.addPropertyChangeListener
            (
                new MenuItemActionChangedListener(menuItem)
            );
            menuItem.setEnabled(action.isEnabled());
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    void clearActions()
    {
        popupMenu_      = null;
        defaultAction_  = null;
        defaultContext_ = null;
    }
//
//    //--------------------------------------------------------------------------
//
//    Action getAction(String actionName)
//    {
//        if (popupMenu_ != null)
//        {
//            MenuElement[] menuElements = popupMenu_.getSubElements();
//            for (int i = 0; i < menuElements.length; i++)
//            {
//                JMenuItem menuItem = menuElements[i].getComponent();
//            }
//        }
//    }
//

    //--------------------------------------------------------------------------

    void invokeDefaultAction()
    {
        try
        {
            if (defaultAction_ != null)
            {
                ActionEvent evt = new ActionEvent
                (
                    this,
                    nodeID_,
                    defaultContext_
                );
                defaultAction_.actionPerformed(evt);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    void showPopupMenu(JComponent parent, int x, int y)
    {
        // Show the popup menu if we have one.
        if (popupMenu_ != null)
        {
            popupMenu_.show(parent, x, y);
        }
    }

    //--------------------------------------------------------------------------
    //----- FoamXTreeRenderer.FoamXTreeItem Interface
    //--------------------------------------------------------------------------

    public String getText()
    {
        return nodeLabel_;
    }
    public Icon getIcon()
    {
        return icon_;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}


