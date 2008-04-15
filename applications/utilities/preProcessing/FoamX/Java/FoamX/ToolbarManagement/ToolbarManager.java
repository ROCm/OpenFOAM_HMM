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
package FoamX.ToolbarManagement;

import java.util.*;

import FoamX.App;
import FoamX.Exceptions.FoamXException;

public class ToolbarManager
    extends javax.swing.JPanel
{
    //--------------------------------------------------------------------------

    /** Horizontal gap between toolbars. */
    static int HGAP = 4;
    /** Vertical gap between toolbars. */
    static int VGAP = 6;

    // Map of toolbar objects, keyed by toolbar name.
    protected java.util.Hashtable toolbars_; 

    //--------------------------------------------------------------------------
    /** Creates new form ToolbarManager */
    public ToolbarManager()
    {
        initComponents();

        toolbars_ = new Hashtable();
    }

    //--------------------------------------------------------------------------

    public boolean isToolbarDefined(String name)
    {
        return toolbars_.containsKey(name);
    }

    //--------------------------------------------------------------------------

    public FoamX.ToolbarManagement.Toolbar getToolbar
    (
        String name
    ) throws FoamXException
    {
        if (!toolbars_.containsKey(name))
        {
            throw new FoamXException("Invalid toolbar name : " + name);
        }

        Toolbar toolbar = (Toolbar)toolbars_.get(name);
        return toolbar;
    }

    //--------------------------------------------------------------------------

    /* create new toolBar */
    public FoamX.ToolbarManagement.Toolbar addToolbar
    (
        String name,
        String description
    ) throws FoamXException
    {
        if (toolbars_.containsKey(name))
        {
            throw new FoamXException("Duplicate toolbar name : " + name);
        }

        // Create new Toolbar object.
        Toolbar toolbar = new Toolbar(name, description, true);
        toolbar.setFont(new java.awt.Font("Dialog", 0, 10));

        // Add to map and content pane.
        toolbars_.put(name, toolbar);
        add(toolbar);
        toolbar.addRef();

        // Refresh the gui.
        updateUI();
        //toolbar.updateUI();

        return toolbar;
    }

    //--------------------------------------------------------------------------

    public void removeToolbar(String name) throws FoamXException
    {
        if (!toolbars_.containsKey(name))
        {
            throw new FoamXException("Invalid toolbar name : " + name);
        }

        Toolbar toolbar = (Toolbar)toolbars_.get(name);
        // Decrement reference count;
        if (toolbar.release() == 0)
        {
            // Remove toolbar from container.
            remove(toolbar);
            updateUI();
        }
    }

    //--------------------------------------------------------------------------

    /* create or reuse existing toolBar */
    public FoamX.ToolbarManagement.Toolbar addMultiToolbar
    (
        String name,
        String description,
        String clientKey
    ) throws FoamXException
    {
        Toolbar toolbar = null;

        if (toolbars_.containsKey(name))
        {
            // If toolbar already exists, increment it's reference count.
            toolbar = (Toolbar)toolbars_.get(name);
            toolbar.addRef();
        }
        else
        {
            // Create new Toolbar object.
            toolbar = new Toolbar(name, description, true);
            toolbar.setFont(new java.awt.Font("Dialog", 0, 10));

            // Add to map and content pane.
            toolbars_.put(name, toolbar);
            add(toolbar);
            toolbar.addRef();
        }

        // Refresh gui.
        updateUI();
        //toolbar.updateUI();

        return toolbar;
    }

    //--------------------------------------------------------------------------

    public void removeMultiToolbar
    (
        String name,
        String clientKey
    ) throws FoamXException
    {
        try
        {
            if (!toolbars_.containsKey(name))
            {
                throw new FoamXException("Invalid toolbar name : " + name);
            }

            // Get specified toolbar.
            Toolbar toolbar = (Toolbar)toolbars_.get(name);

            // Remove mutil-actions.
            toolbar.removeMultiActions(clientKey);

            // Decrement reference count;
            if (toolbar.release() <= 0)
            {
                // Remove toolbar from container.
                toolbars_.remove(name);
                remove(toolbar);
                updateUI();
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents()//GEN-BEGIN:initComponents
    {
        
        setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 2, 2));
        
        setBorder(new javax.swing.border.EtchedBorder());
    }//GEN-END:initComponents


    //--------------------------------------------------------------------------

    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
