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

import java.net.URL;
import java.util.Hashtable;
import java.awt.Component;
import java.awt.event.ActionEvent;
import javax.swing.*;
import javax.swing.tree.*;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;

import FoamX.App;
import FoamX.Util.FoamXTreeRenderer;

public class PatchPhysicalTypesTab
    extends javax.swing.JPanel
{
    //--------------------------------------------------------------------------

    private BoundaryDefinitionModel boundaryDefinitionModel_;
    private BoundaryDefinitionPanel boundaryDefinitionPanel_;
    private AddPatchPhysicalTypeAction addPatchPhysicalTypeAction_;
    private DeletePatchPhysicalTypeAction deletePatchPhysicalTypeAction_;

    //--------------------------------------------------------------------------
    /** Constructor for DictionaryPanel. */
    public PatchPhysicalTypesTab(ApplicationModel appModel)
    {
        try
        {
            boundaryDefinitionModel_ = appModel.getBoundaryDefinitionModel();
            boundaryDefinitionPanel_ = new BoundaryDefinitionPanel(appModel);

            // Initialise the gui.
            initComponents();

            // Set tree renderer.
            tree_.setCellRenderer(new FoamXTreeRenderer());
            tree_.putClientProperty("JTree.lineStyle", "Angled");

            // Set single selection mode.
            tree_.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);

            // Listen out for tree node selection events.
            tree_.getSelectionModel().addTreeSelectionListener
            (
                new TreeSelectionListener()
                {
                    public void valueChanged(TreeSelectionEvent evt)
                    {
                        OnSelectionChanged(evt);
                    }
                }
            );

            // Make sure the tree is displayed properly.
            tree_.updateUI();

            // Initialise the toolbar. Create a button for each action.
            addPatchPhysicalTypeAction_    = new AddPatchPhysicalTypeAction();
            deletePatchPhysicalTypeAction_ = new DeletePatchPhysicalTypeAction();
            addToolbarButton(addPatchPhysicalTypeAction_);
            addToolbarButton(deletePatchPhysicalTypeAction_);

            // Initialise actions
            addPatchPhysicalTypeAction_.setEnabled(false);
            deletePatchPhysicalTypeAction_.setEnabled(false);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    private void addToolbarButton(Action action)
    {
        // Add new button to the toolbar using the action object.
        javax.swing.JButton button = toolbar_.add(action);

        // Set the tooltip text.
        button.setToolTipText((String)action.getValue(Action.SHORT_DESCRIPTION));
        button.setText("");
        button.setFont(toolbar_.getFont());      // Use same font as toolbar.
    }

    //--------------------------------------------------------------------------

    public void updateModel()
    {
    }

    //--------------------------------------------------------------------------

    private String getNewBoundaryDefinitionName()
    {
        // Prompt user for a new Boundary Definition name.
        return (String)JOptionPane.showInputDialog(this,
                                                    "Enter New Boundary Definition Name",
                                                    "Boundary Definition Name",
                                                    JOptionPane.PLAIN_MESSAGE, null, null, null);
    }

    //--------------------------------------------------------------------------
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        contextMenu_ = new javax.swing.JPopupMenu();
        subPanel_ = new javax.swing.JPanel();
        splitPane_ = new javax.swing.JSplitPane();
        scrollPaneRight_ = new javax.swing.JScrollPane();
        leftPanel_ = new javax.swing.JPanel();
        buttonPanel_ = new javax.swing.JPanel();
        toolbar_ = new javax.swing.JToolBar();
        scrollPaneTree_ = new javax.swing.JScrollPane();
        tree_ = new javax.swing.JTree();
        
        contextMenu_.setFont(new java.awt.Font("Dialog", 0, 10));
        
        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        subPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        splitPane_.setDividerLocation(200);
        splitPane_.setLastDividerLocation(200);
        scrollPaneRight_.setAutoscrolls(true);
        scrollPaneRight_.setViewportView( boundaryDefinitionPanel_ );
        
        splitPane_.setRightComponent(scrollPaneRight_);
        
        leftPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints3;
        
        buttonPanel_.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT));
        
        buttonPanel_.setBorder(new javax.swing.border.EtchedBorder());
        toolbar_.setFloatable(false);
        toolbar_.setFont(new java.awt.Font("Dialog", 0, 8));
        buttonPanel_.add(toolbar_);
        
        gridBagConstraints3 = new java.awt.GridBagConstraints();
        gridBagConstraints3.gridx = 0;
        gridBagConstraints3.gridy = 0;
        gridBagConstraints3.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints3.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints3.weightx = 1.0;
        leftPanel_.add(buttonPanel_, gridBagConstraints3);
        
        scrollPaneTree_.setAutoscrolls(true);
        tree_.setFont(new java.awt.Font("Dialog", 0, 10));
        tree_.setShowsRootHandles(true);
        tree_.setModel( boundaryDefinitionModel_.getTreeModel() );
        tree_.setAutoscrolls(true);
        tree_.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                OnTreeMouseClicked(evt);
            }
        });
        
        scrollPaneTree_.setViewportView(tree_);
        
        gridBagConstraints3 = new java.awt.GridBagConstraints();
        gridBagConstraints3.gridx = 0;
        gridBagConstraints3.gridy = 1;
        gridBagConstraints3.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints3.anchor = java.awt.GridBagConstraints.SOUTH;
        gridBagConstraints3.weightx = 1.0;
        gridBagConstraints3.weighty = 1.0;
        leftPanel_.add(scrollPaneTree_, gridBagConstraints3);
        
        splitPane_.setLeftComponent(leftPanel_);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints2.weightx = 1.0;
        gridBagConstraints2.weighty = 1.0;
        subPanel_.add(splitPane_, gridBagConstraints2);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.insets = new java.awt.Insets(10, 10, 10, 10);
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        add(subPanel_, gridBagConstraints1);
        
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

  private void OnTreeMouseClicked (java.awt.event.MouseEvent evt)
    {//GEN-FIRST:event_OnTreeMouseClicked
        try
        {
            // Get the node object for the current selection.
            TreePath tp = tree_.getPathForLocation(evt.getX(), evt.getY());
            if (tp != null)
            {
                // Check for popup menu trigger.
                //if (evt.isPopupTrigger())
                if (evt.getModifiers() == evt.BUTTON3_MASK)   // Check explicitly for right mouse click.
                {
                    DefaultMutableTreeNode      nodeItem = (DefaultMutableTreeNode)tp.getLastPathComponent();
                    BoundaryDefinitionModelItem nodeInfo = (BoundaryDefinitionModelItem)nodeItem.getUserObject();

                    // Construct an appropriate context menu,
                    contextMenu_.removeAll();
                    java.awt.Font font = contextMenu_.getFont();

                    if (nodeInfo.isPatchType())
                    {
                        contextMenu_.add(addPatchPhysicalTypeAction_).setFont(font);
                    }
                    else if (nodeInfo.isBoundaryDefinition())
                    {
                        contextMenu_.add(addPatchPhysicalTypeAction_).setFont(font);
                        contextMenu_.add(deletePatchPhysicalTypeAction_).setFont(font);
                    }

                    // Show context menu.
                    contextMenu_.show(tree_, evt.getX(), evt.getY());
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
  }//GEN-LAST:event_OnTreeMouseClicked

    //--------------------------------------------------------------------------

    public void OnSelectionChanged(TreeSelectionEvent evt)
    {
        try
        {
            // Get the node object for the current selection.
            if (evt.isAddedPath())    // Has the new node been added to the selection?
            {
                TreePath tp = evt.getPath();
                if (tp != null)
                {
                    // Select this node in the model.
                    DefaultMutableTreeNode      nodeItem = (DefaultMutableTreeNode)tp.getLastPathComponent();
                    BoundaryDefinitionModelItem nodeInfo = (BoundaryDefinitionModelItem)nodeItem.getUserObject();

                    // Set the action enabled/disabled status.
                    addPatchPhysicalTypeAction_.setEnabled(nodeInfo.isPatchType() || nodeInfo.isBoundaryDefinition());
                    deletePatchPhysicalTypeAction_.setEnabled(nodeInfo.isBoundaryDefinition() && nodeItem.getChildCount() == 0); // Can't delete anything with a child, yet.

                    boundaryDefinitionModel_.select(nodeItem);
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPopupMenu contextMenu_;
    private javax.swing.JPanel subPanel_;
    private javax.swing.JSplitPane splitPane_;
    private javax.swing.JScrollPane scrollPaneRight_;
    private javax.swing.JPanel leftPanel_;
    private javax.swing.JPanel buttonPanel_;
    private javax.swing.JToolBar toolbar_;
    private javax.swing.JScrollPane scrollPaneTree_;
    private javax.swing.JTree tree_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------

    class AddPatchPhysicalTypeAction
    extends AbstractAction
    {
        AddPatchPhysicalTypeAction()
        {
            putValue(Action.SMALL_ICON, App.getResources().getIcon("NewPatchPhysicalTypeImage"));
            putValue(Action.NAME, "Add Boundary Type");
            putValue(Action.SHORT_DESCRIPTION, "Add Boundary Type");
            putValue(Action.LONG_DESCRIPTION, "Add Boundary Type");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Add a new boundary type to the current selection.
            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // See if this is a patch node or a boundary definition node.
                    DefaultMutableTreeNode      parentNode = (DefaultMutableTreeNode)tp.getLastPathComponent();
                    BoundaryDefinitionModelItem nodeInfo   = (BoundaryDefinitionModelItem)parentNode.getUserObject();
                    // Make sure the selected node is not the root node.
                    if (!nodeInfo.isRoot())
                    {
                        String name = getNewBoundaryDefinitionName();
                        if (name != null && name.length()> 0)
                        {
                            DefaultMutableTreeNode newNode = boundaryDefinitionModel_.addBoundaryDefinitionNode(name, parentNode);
                            if (newNode != null)
                            {
                                // If addition successful, select the new child node.
                                tree_.expandPath(tp);
                                TreePath childTP = new TreePath(newNode.getPath());
                                tree_.getSelectionModel().clearSelection();
                                tree_.getSelectionModel().setSelectionPath(childTP);
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    class DeletePatchPhysicalTypeAction
    extends AbstractAction
    {
        DeletePatchPhysicalTypeAction()
        {
            putValue(Action.SMALL_ICON, App.getResources().getIcon("DelPatchPhysicalTypeImage"));
            putValue(Action.NAME, "Delete Boundary Type");
            putValue(Action.SHORT_DESCRIPTION, "Delete Boundary Type");
            putValue(Action.LONG_DESCRIPTION, "Delete Boundary Type");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Delete the selected boundary definition.
            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // See if this is a patch node or a boundary definition node with children.
                    DefaultMutableTreeNode      nodeItem = (DefaultMutableTreeNode)tp.getLastPathComponent();
                    BoundaryDefinitionModelItem nodeInfo = (BoundaryDefinitionModelItem)nodeItem.getUserObject();
                    // Make sure the selected node is a boundary type node.
                    if (nodeInfo.isBoundaryDefinition())
                    {
                        if (boundaryDefinitionModel_.deleteBoundaryDefinitionNode(nodeItem))
                        {
                            // If deletion successful, select the parent node.
                            tree_.getSelectionModel().clearSelection();
                            tree_.getSelectionModel().setSelectionPath(tp.getParentPath());
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}



