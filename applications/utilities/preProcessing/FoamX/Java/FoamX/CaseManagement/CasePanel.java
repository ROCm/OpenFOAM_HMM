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

import java.io.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.FocusListener;
import java.awt.event.FocusEvent;
import java.util.*;
import javax.swing.*;
import javax.swing.tree.*;
import javax.swing.event.*;

import FoamX.App;
import FoamX.Exceptions.FoamXException;
import FoamX.CaseManagement.HostChooserDlg;
import FoamX.Reporting.ReportingWindow;
import FoamX.ToolbarManagement.Toolbar;
import FoamX.Tools.InvokeUtilityAction;
import FoamX.Util.BusyCursor;
import FoamX.Util.CollapseTreeAction;
import FoamX.Util.ExpandTreeAction;
import FoamX.Util.FoamXTreeRenderer;
import FoamX.Util.TabSelection;
import FoamX.Util.RunPanel;
import FoamX.Util.FileEvent;
import FoamX.WindowManagement.FoamXInternalFrame;


import FoamXServer.CaseServer.ICaseServer;
import FoamXServer.CaseBrowser.ICaseBrowser;
import FoamXServer.IDictionaryEntry;
import FoamXServer.IDictionaryEntryHolder;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;
import FoamXServer.ApplicationDescriptor;

public class CasePanel
    extends javax.swing.JPanel
    implements IModuleHost, TabSelection
{
    //--------------------------------------------------------------------------

    public static final boolean KILL_SERVER = true;
    public static final boolean KEEP_SERVER = false;


    protected static final String TOOLBAR_NAME  = "CaseToolbar";

    protected ICaseServer caseServer_;
    protected ICaseBrowser caseBrowser_;
    protected EventListenerList listenerList_;
    protected Hashtable moduleMap_;       // Map of IModule objects.
    protected Hashtable caseActionMap_;   // Map of possible actions
    protected DefaultTreeModel treeModel_;
    protected DefaultMutableTreeNode rootNode_;
    protected Toolbar toolBar_;
    protected String toolBarKey_;
    protected int moduleCount_;
    protected int treeNodeCount_;

    // Case Panel actions.
    CloseCaseAction closeCaseAction_;
    SaveCaseAction saveCaseAction_;
    StartCalculationAction startCalculationAction_;
    RunCalculationAction runCalculationAction_;
    StopCalculationAction stopCalculationAction_;
    ExpandTreeAction expandAllAction_;
    CollapseTreeAction collapseAllAction_;
    InvokeUtilityAction invokeUtilityAction_;

    //--------------------------------------------------------------------------
    /** CasePanel constructor. */
    public CasePanel(ICaseBrowser caseBrowser, ICaseServer caseServer)
        throws Exception
    {
        try
        {
            // Initialise the GUI components.
            initComponents();

            // Create listener list.
            listenerList_ = new EventListenerList();

            // Create module map.
            moduleMap_ = new Hashtable();

            // Create case action map
            caseActionMap_ = new Hashtable();

            // Store the primary reference to the case server.
            caseServer_ = caseServer;
            caseServer_.managed(true);   // Mark this case as managed.

            // Store the primary reference to the case browser.
            // (is used for invoking external actions e.g. running)
            caseBrowser_ = caseBrowser;

            // Initialise the tree model.
            Icon icon = App.getResources().getIcon("Case.CaseRootImage");
            TreeNodeInfo nodeInfo =
                new TreeNodeInfo(caseServer_.caseName(), 0, icon);
            rootNode_      = new DefaultMutableTreeNode(nodeInfo, true);
            treeModel_     = new DefaultTreeModel(rootNode_);
            treeNodeCount_ = 1;

            // Set the tree model and renderer.
            tree_.setModel(treeModel_);
            tree_.setCellRenderer(new FoamXTreeRenderer());
            tree_.putClientProperty
            ( 
                "JTree.lineStyle",
                App.getResources().getResourceString("Tree.LineStyle")
            );

            // Initialise actions.
            closeCaseAction_        = new CloseCaseAction();
            saveCaseAction_         = new SaveCaseAction();
            runCalculationAction_   = new RunCalculationAction();
            startCalculationAction_ = new StartCalculationAction();
            stopCalculationAction_  = new StopCalculationAction();
            expandAllAction_        = new ExpandTreeAction(tree_);
            collapseAllAction_      = new CollapseTreeAction(tree_);
            invokeUtilityAction_ =
                new InvokeUtilityAction
                (
                    caseBrowser_,
                    caseServer_.caseRoot(),
                    caseServer_.caseName()
                );

            // Initialise toolbar. Toolbar may already exist.
            toolBarKey_ =
                caseServer_.caseRoot() + '/' + caseServer_.caseName();
            toolBar_ =
                App.getToolbarManager().addMultiToolbar
                (
                    TOOLBAR_NAME,
                    "Case Toolbar",
                    toolBarKey_
                );

            addToolbarButton(closeCaseAction_);
            addToolbarButton(saveCaseAction_);
            toolBar_.addSeparator();
            addToolbarButton(runCalculationAction_);
            addToolbarButton(startCalculationAction_);
            addToolbarButton(stopCalculationAction_);
            toolBar_.addSeparator();
            addToolbarButton(expandAllAction_);
            addToolbarButton(collapseAllAction_);

            // Attach utilities context menu to root node.
            buildUtilsMenu();

            // Reset all actions.
            resetActions(true);

            // Initialise and activate the modules for this case.
            loadCaseModules();
            activateCaseModules();

            // Listen out for selection events.
            tree_.getSelectionModel().addTreeSelectionListener
            (
                new TreeSelectionListener()
                {
                    public void valueChanged(TreeSelectionEvent evt)
                    {
                        OnTreeSelectionChanged(evt);
                    }
                }
            );

            // Expand the root node.
            TreePath tp = new TreePath(treeModel_.getPathToRoot(rootNode_));
            if (tp != null)
            {
                tree_.expandPath(tp);
                tree_.updateUI();
            }

            // Print success message to FoamX tab in reporting window
            printAppMessage
            (
                "Opened case " + caseServer_.caseName() + 
                " at database time " + caseServer_.getTime()
            );

            // Cannot print to case tab in reporting window since is not
            // yet open!
        }
        catch (FoamXIOError ioErr)
        {
            // Problem: caseManager hasn't yet registered me so cannot use
            // closeCaseAction_.actionPerformed(null)
            // So close myself and throw up to caseManager.
            caseServer_.managed(false);
            closeCase(KEEP_SERVER);
            // Handle at higher levels
            throw ioErr;
        }
        catch (Exception ex)
        {
            App.handleException(ex);
        }
    }

    //--------------------------------------------------------------------------

    public String getCaseRoot()
    {
        return caseServer_ == null ? "" : caseServer_.caseRoot();
    }

    //--------------------------------------------------------------------------

    public String getCaseName()
    {
        return caseServer_ == null ? "" : caseServer_.caseName();
    }

    //--------------------------------------------------------------------------
    /** Close the case and release the case server object.
     *  Argument determines whether the caseServer gets killed. Is tricky.
     *  Is usually true, except when CasePanel fails to construct
     *  correctly and caseServer.managed(false) is used to signal
     *  the caseBrowserPanel that it failed.
     *  (the normal cleanup action through the
     *  status listeners does not function if CasePanel has not been
     *  constructed correctly)
     */
    void closeCase(boolean killServer)
    {
        try
        {
            // Shutdown all case modules.
            shutdownCaseModules();

            // Close down the case server and unlock the case.
            if (killServer && (caseServer_ != null))
            {
                caseServer_.close();
                caseServer_ = null;
            }

            // Remove toolbar.
            App.getToolbarManager().removeMultiToolbar
            (
                TOOLBAR_NAME,
                toolBarKey_
            );
            toolBar_    = null;
            toolBarKey_ = "";
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }


    //--------------------------------------------------------------------------
    /** Close the case but ask whether to save case beforehand */
    void closeCaseNice(boolean killServer)
    {
        try
        {
            // Check if should have been saved.
            if (caseServer_ != null)
            {
                checkAndSave
                (
                    false,
                    caseServer_.caseName() + " : case has been changed."
                    + " Save before closing?",
                    "Save Case"
                );
            }

            // Close case as normal
            closeCase(killServer);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /**
      * helper function to save case dialog if needed
      * @return false if cancel pressed, true otherwise
      */
    protected boolean checkAndSave
    (
        boolean allowCancel,
        String msg,
        String title
    )
    {
        boolean cancelled = false;

        if (caseServer_.modified())
        {
            int option;
            if (allowCancel)
            {
                option = JOptionPane.YES_NO_CANCEL_OPTION;
            }
            else
            {
                option = JOptionPane.YES_NO_OPTION;
            }

            // Prompt the user to see if the case files are
            // to be saved.
            int ret = JOptionPane.showConfirmDialog
            (
                App.getRootFrame(),
                msg,
                title,
                option,
                JOptionPane.QUESTION_MESSAGE
            );
            if (ret == JOptionPane.OK_OPTION)
            {
                // Behave as if 'save' button pressed
                saveCaseAction_.actionPerformed(null);
            }
            else if (ret == JOptionPane.CANCEL_OPTION)
            {
                cancelled = true;
            }
        }
        return !cancelled;
    }

    //--------------------------------------------------------------------------

    protected void printAppMessage(String message)
    {
        // Get reporting window for Application.
        ReportingWindow reportWin =
            App.getReportingManager().getReportingWindow("FoamX");
        if (reportWin != null)
        {
            // Print message into main reporting window.
            reportWin.printMessage(message);
        }
    }

    //--------------------------------------------------------------------------

    protected void printCaseMessage(String message)
    {
        // Get reporting window for this case.
        String reportWinName = getCaseRoot() + '/' + getCaseName();
        ReportingWindow reportWin =
            App.getReportingManager().getReportingWindow(reportWinName);
        if (reportWin != null)
        {
            // Print message into main reporting window.
            reportWin.printMessage(message);
        }
    }

    //--------------------------------------------------------------------------

    protected DefaultMutableTreeNode getNode(int nodeID)
    {
        DefaultMutableTreeNode node = null;

        try
        {
            // Find the specified node.
            Enumeration nodeIter = rootNode_.depthFirstEnumeration();
            while (nodeIter.hasMoreElements())
            {
                DefaultMutableTreeNode nodeItem =
                    (DefaultMutableTreeNode)nodeIter.nextElement();
                TreeNodeInfo fxNode = (TreeNodeInfo)nodeItem.getUserObject();
                if (fxNode.getNodeID() == nodeID)
                {
                    node = nodeItem;
                    break;
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return node;
    }

    //--------------------------------------------------------------------------

    protected TreeNodeInfo getNodeInfo(int nodeID)
    {
        TreeNodeInfo fxNode = null;

        try
        {
            // Find the specified node.
            DefaultMutableTreeNode node = getNode(nodeID);
            if (node != null)
            {
                fxNode = (TreeNodeInfo)node.getUserObject();
            }
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
            App.handleAllExceptions(ex);
        }

        return fxNode;
    }

    //--------------------------------------------------------------------------

    protected void loadCaseModules() throws FoamXIOError
    {
        try
        {
            // Get the list of modules to load from the application class
            // object.
            String[] moduleNames = caseServer_.application().modules();

            // Load and initialise all modules.
            for (short i = 0; i <moduleNames.length; i++)
            {
                // Only initialise this module if we have not added it before.
                if (!moduleMap_.containsKey(moduleNames[i]))
                {
                    IModule module = initialiseModule(moduleNames[i]);
                    if (module != null)
                    {
                        moduleMap_.put(moduleNames[i], module);
                    }
                }
            }
        }
        catch (FoamXIOError ioErr)
        {
            // C++ IO error: handle at higher levels
            throw ioErr;
        }
        catch (Exception ex)
        {
            App.handleException(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Create a FoamXModule and initialise it. */
    protected IModule initialiseModule(String moduleName)
       throws FoamXIOError
    {
        IModule module = null;

        try
        {
            // Load the specified module.
            ClassLoader loader = ClassLoader.getSystemClassLoader();
            Object obj = java.beans.Beans.instantiate(loader, moduleName);

            // Make sure that this module supports the IModule interface.
            java.lang.Class c;
            c = Class.forName("FoamX.CaseManagement.IModule");
            if (!java.beans.Beans.isInstanceOf(obj, c))
            {
                throw new Exception
                (
                    "Module " + moduleName
                    + " does not support the IModule interface."
                );
            }

            // Get the IModule interface.
            module = (IModule)obj;

            // Allocate a new module ID.
            moduleCount_++;

            // Finally, call the module's Initialise method.
            if (!module.initialiseModule(this, moduleCount_))
            {
                throw new Exception
                (
                    "Module initialisation failure for module " + moduleName
                    + "."
                );
            }
        }
        catch (FoamXIOError ioErr)
        {
            // C++ IO error: handle at higher levels
            throw ioErr;
        }
        catch (FoamXError fxErr)
        {
            // Corba error: handle here. What to do with these?
            App.handleException(fxErr);
        }
        catch (Exception ex)
        {
            // Other: handle here
            App.handleException(ex);
        }

        return module;
    }

    //--------------------------------------------------------------------------

    protected void shutdownCaseModules()
    {
        try
        {
            // Shutdown all case modules.
            Enumeration iter = moduleMap_.elements();
            while (iter.hasMoreElements())
            {
                // Shutdown the module.
                IModule module = (IModule)iter.nextElement();
                module.shutdownModule();
            }
            // Clear the module map.
            moduleMap_.clear();
            // Clear the actions map
            caseActionMap_.clear();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    protected void activateCaseModules()
    {
        try
        {
            // Activate all case modules.
            Enumeration iter = moduleMap_.elements();
            while (iter.hasMoreElements())
            {
                // Shutdown the module.
                IModule module = (IModule)iter.nextElement();
                module.activate();
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    protected void deactivateCaseModules()
    {
        try
        {
            // Activate all case modules.
            Enumeration iter = moduleMap_.elements();
            while (iter.hasMoreElements())
            {
                // Shutdown the module.
                IModule module = (IModule)iter.nextElement();
                module.deactivate();
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Reset all actions to the default (disabled) state. */
    protected void resetActions(boolean tabSelected)
    {
        closeCaseAction_.setEnabled(false);
        saveCaseAction_.setEnabled(false);
        runCalculationAction_.setEnabled(false);
        startCalculationAction_.setEnabled(false);
        stopCalculationAction_.setEnabled(false);
        expandAllAction_.setEnabled(false);
        collapseAllAction_.setEnabled(false);
        invokeUtilityAction_.setEnabled(false);
    }

    //--------------------------------------------------------------------------
    /** Set all actions according to the current selection. */
    protected void setActions()
    {
        try
        {
            resetActions(true);

            closeCaseAction_.setEnabled(true);
            saveCaseAction_.setEnabled(true);
            runCalculationAction_.setEnabled(true);
            startCalculationAction_.setEnabled(true);
            stopCalculationAction_.setEnabled(true);
            expandAllAction_.setEnabled(true);
            collapseAllAction_.setEnabled(true);
            invokeUtilityAction_.setEnabled(true);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    protected void addToolbarButton(AbstractAction action)
    {
        // Add new button to the toolbar using the action object.
        String name = (String)action.getValue(Action.NAME);
        javax.swing.JButton button =
            toolBar_.addMultiAction(name, toolBarKey_, action);

        // Set the tooltip text.
        button.setToolTipText
        (
            (String)action.getValue(Action.SHORT_DESCRIPTION)
        );
        button.setText("");
        button.setFont(toolBar_.getFont());      // Use same font as toolbar.
    }

    //--------------------------------------------------------------------------

    protected void buildUtilsMenu()
    {
        try
        {
            Hashtable utilMenuMap = new Hashtable(5);

            // Construct utilities context menu.
            JPopupMenu contextMenu = new JPopupMenu();
            contextMenu.setFont(new java.awt.Font("Dialog", 0, 10));

            JMenu utilsMenu = new JMenu("Foam Utilities");
            utilsMenu.setFont(contextMenu.getFont());
            contextMenu.add(utilsMenu);

            // Get utility descriptors.
            ApplicationDescriptor[] utilities =
                caseServer_.foamProperties().utilities();

            for (int i = 0; i < utilities.length; i++)
            {
                // Build category sub-menus.
                JMenu menu = utilsMenu;
                if
                (
                    (utilities[i].category != null)
                 && (utilities[i].category.length() > 0)
                )
                {
                    // Add all required sub-menus.
                    StringTokenizer tok =
                        new StringTokenizer(utilities[i].category, "/");
                    String catKey = new String();
                    while (tok.hasMoreTokens())
                    {
                        String category = tok.nextToken();
                        catKey = catKey + category;
                        if (utilMenuMap.containsKey(catKey))
                        {
                            menu = (JMenu)utilMenuMap.get(catKey);
                        }
                        else
                        {
                            JMenu catMenu = new JMenu(category);
                            catMenu.setFont(contextMenu.getFont());
                            menu.add(catMenu);
                            utilMenuMap.put(catKey, catMenu);
                            menu = catMenu;
                        }
                    }
                }

                // Add a menu item for this utility.
                JMenuItem utilItem = new JMenuItem(utilities[i].name);
                utilItem.setFont(contextMenu.getFont());

                // Hook up to invoke utility action.
                utilItem.addActionListener(invokeUtilityAction_);

                // Specify which utility to invoke.
                utilItem.putClientProperty
                (
                    "utilityDescriptor",
                    utilities[i]
                );

                // Add to menu.
                menu.add(utilItem);
            }

            // Attach context menu to root node.
            TreeNodeInfo nodeInfo = (TreeNodeInfo)rootNode_.getUserObject();
            nodeInfo.setContextMenu(contextMenu);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //---- TabSelection Interface Methods
    //--------------------------------------------------------------------------

    public void tabSelected()
    {
        // Check to see if we still have a toolbar.
        if (toolBar_ != null)
        {
            // Get the muti-action toolbar to direct the actions to this object.
            toolBar_.selectClient(toolBarKey_);
            setActions();
        }
    }

    public void tabDeselected()
    {
        resetActions(false);
    }

    //--------------------------------------------------------------------------
    //---- CaseStatusListener Methods
    //--------------------------------------------------------------------------

    public void addCaseStatusListener(CaseStatusListener l)
    {
        listenerList_.add(CaseStatusListener.class, l);
    }

    //--------------------------------------------------------------------------

    public void removeCaseStatusListener(CaseStatusListener l)
    {
        listenerList_.remove(CaseStatusListener.class, l);
    }

    //--------------------------------------------------------------------------

    protected void fireCloseCase(String caseRoot, String caseName)
    {
        // Create event object.
        CaseStatusEvent evt =
            new CaseStatusEvent(this, caseRoot, caseName, caseServer_);

        // Process the listeners last to first, notifying those that
        // are interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseStatusListener.class)
            {
                ((CaseStatusListener)listeners[i+1]).caseClosed(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents()//GEN-BEGIN:initComponents
    {
        contextMenu_ = new javax.swing.JPopupMenu();
        scrollPanel_ = new javax.swing.JScrollPane();
        tree_ = new javax.swing.JTree();
        
        
        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        tree_.setFont(new java.awt.Font("Dialog", 0, 10));
        tree_.setShowsRootHandles(true);
        tree_.setToggleClickCount(3);
        tree_.addMouseListener(new java.awt.event.MouseAdapter()
        {
            public void mouseClicked(java.awt.event.MouseEvent evt)
            {
                OnMouseClicked(evt);
            }
        });
        
        scrollPanel_.setViewportView(tree_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.SOUTHWEST;
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        add(scrollPanel_, gridBagConstraints1);
        
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

  private void OnMouseClicked (java.awt.event.MouseEvent evt)
    {//GEN-FIRST:event_OnMouseClicked
        try
        {
            // Get the node object for the current selection.
            TreePath tp = tree_.getPathForLocation(evt.getX(), evt.getY());
            if (tp != null)
            {
                // Make sure the node is selected.
                tree_.getSelectionModel().clearSelection();
                tree_.getSelectionModel().setSelectionPath(tp);

                DefaultMutableTreeNode nodeItem =
                    (DefaultMutableTreeNode)tp.getLastPathComponent();
                TreeNodeInfo fxNode =
                    (TreeNodeInfo)nodeItem.getUserObject();

                // Check for popup menu trigger.
                //else if (evt.isPopupTrigger())
                if (evt.getModifiers() == evt.BUTTON3_MASK)
                {
                    if (fxNode.getContextMenu() != null)
                    {
                        fxNode.showPopupMenu
                        (
                            (JComponent)evt.getSource(),
                            evt.getX(),
                            evt.getY()
                        );
                    }
                }
                else if (evt.getClickCount() == 2)
                {
                    // If there is a default command bound to this node,
                    // invoke it.
                    fxNode.invokeDefaultAction();
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
  }//GEN-LAST:event_OnMouseClicked

    //--------------------------------------------------------------------------

    public void OnTreeSelectionChanged(TreeSelectionEvent evt)
    {
        try
        {
            // Get the node object for the current selection.

            // Has the new node been added to the selection?
            if (evt.isAddedPath())
            {
                TreePath tp = evt.getPath();
                if (tp != null)
                {
                    setActions();
                    //DefaultMutableTreeNode nodeItem = (DefaultMutableTreeNode)tp.getLastPathComponent();
                    //TreeNodeInfo           fxNode   = (TreeNodeInfo)nodeItem.getUserObject();
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
    private javax.swing.JScrollPane scrollPanel_;
    private javax.swing.JTree tree_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //---- IModuleHost Interface
    //--------------------------------------------------------------------------

    public JFrame getRootFrame()
    {
        // Return root frame object.
        return App.getRootFrame();
    }

    //--------------------------------------------------------------------------

    public ICaseServer getCaseServer()
    {
        // Return our case server reference.
        return caseServer_;
    }

    //--------------------------------------------------------------------------

    public ICaseBrowser getCaseBrowser()
    {
        return caseBrowser_;
    }

    //--------------------------------------------------------------------------

    public void printMessage(String message)
    {
        // Get the window key for this case.
        String windowKey =
            caseServer_.caseRoot() + '/' + caseServer_.caseName();

        // Get the ReportingWindow object for this case.
        ReportingWindow reportingWindow =
            App.getReportingManager().getReportingWindow(windowKey);

        if (reportingWindow != null)
        {
            reportingWindow.printMessage(message);
        }
    }

    //--------------------------------------------------------------------------

    public void setRootLabel(String label)
    {
        // Set the label of the root node.
        TreeNodeInfo fxNode = (TreeNodeInfo)rootNode_.getUserObject();
        fxNode.setNodeLabel(label);
        treeModel_.nodeChanged(rootNode_);
    }

    //--------------------------------------------------------------------------

    public String getRootLabel()
    {
        // Set the label of the root node.
        TreeNodeInfo fxNode = (TreeNodeInfo)rootNode_.getUserObject();
        return fxNode.getNodeLabel();
    }

    //--------------------------------------------------------------------------

    public int addNode(int parentNodeID, String nodeLabel, Icon icon)
    {
        int nodeID = -1;

        try
        {
            if (icon == null)
            {
                throw new FoamXException("Null icon in addNode.");
            }

            // Find the specified parent node.
            DefaultMutableTreeNode parentNode = getNode(parentNodeID);
            if (parentNode == null)
            {
                throw new FoamXException("Invalid parent node ID.");
            }

            // Allocate a new node ID.
            nodeID = treeNodeCount_++;

            // Create new node.
            DefaultMutableTreeNode newNode =
                new DefaultMutableTreeNode
                (
                    new TreeNodeInfo(nodeLabel, nodeID, icon)
                );

            // Insert node into tree.
            treeModel_.insertNodeInto
            (
                newNode,
                parentNode,
                parentNode.getChildCount()
            );
            treeModel_.nodeChanged(parentNode);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
            nodeID = -1;
        }

        return nodeID;
    }

    //--------------------------------------------------------------------------

    public void removeNode(int nodeID)
    {
        try
        {
            // Find the specified node info object.
            TreeNodeInfo fxNode = getNodeInfo(nodeID);
            if (fxNode == null)
            {
                throw new Exception("Invalide node ID.");
            }
            // And clear it.
            fxNode.clearActions();

            // Find the specified node.
            DefaultMutableTreeNode node = getNode(nodeID);
            if (node == null)
            {
                throw new FoamXException("Invalid node ID.");
            }

            // Remove node from tree.
            treeModel_.removeNodeFromParent(node);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void addNodeAction(int nodeID, Action action, String context)
    {
        try
        {
            // Find the specified node info object.
            TreeNodeInfo fxNode = getNodeInfo(nodeID);
            if (fxNode == null)
            {
                throw new Exception("Invalide node ID.");
            }

            // Store reference to action under its name
            caseActionMap_.put(action.getValue(Action.NAME), action);

            // Attach the command to the node.
            fxNode.addAction(action, context);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public Action getAction(String actionName)
    {
        // Return stored action
        return (Action)caseActionMap_.get(actionName);
    }

    //--------------------------------------------------------------------------

    public int createWindow(javax.swing.JComponent component)
    {
        // Create a new window and return it's ID.
        return App.getWindowManager().createWindow(component);
    }

    //--------------------------------------------------------------------------

    public FoamXInternalFrame getWindowFrame(int windowID)
    {
        // Get the window frame from the window manager and return.
        return App.getWindowManager().getWindowFrame(windowID);
    }

    //--------------------------------------------------------------------------

    public void closeWindow(int windowID)
    {
        // Close the specified window.
        App.getWindowManager().closeWindow(windowID);
    }

    //--------------------------------------------------------------------------

    public void showWindow(int windowID)
    {
        // Show the specified window.
        App.getWindowManager().showWindow(windowID);
    }

    //--------------------------------------------------------------------------

    public void hideWindow(int windowID)
    {
        // Hide the specified window.
        App.getWindowManager().hideWindow(windowID);
    }

    //--------------------------------------------------------------------------
    //---- Action Classes
    //--------------------------------------------------------------------------

    private class CloseCaseAction
    extends AbstractAction
    {
        CloseCaseAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("CloseCaseImage")
            );
            putValue(Action.NAME, "Close Case");
            putValue(Action.SHORT_DESCRIPTION, "Close Case");
            putValue(Action.LONG_DESCRIPTION, "Close Case");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(App.getRootFrame());

            try
            {
                if (caseServer_ != null)
                {
                    boolean doClose = checkAndSave
                        (
                            true,
                            caseServer_.caseName() + " : case has been changed."
                            + " Save before closing?",
                            "Save Case"
                        );
                    if (doClose)
                    {
                        // Get the case manager to do all the work.
                        // Case is actually closed in the closeCase method.
                        fireCloseCase
                        (
                            caseServer_.caseRoot(),
                            caseServer_.caseName()
                        );
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

    private class SaveCaseAction
    extends AbstractAction
    {
        SaveCaseAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("SaveCaseImage")
            );
            putValue(Action.NAME, "Save Case");
            putValue(Action.SHORT_DESCRIPTION, "Save Case");
            putValue(Action.LONG_DESCRIPTION, "Save Case");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(App.getRootFrame());

            try
            {
                if (caseServer_ != null)
                {
                    // Temporarily have modules editing this case do nothing
                    deactivateCaseModules();

                    // Save the case.
                    caseServer_.save();

                    // Reactivate modules
                    activateCaseModules();

                    // Print success message.
                    printCaseMessage
                    (
                        "Case " + toolBarKey_ + " saved successfully to time "
                        + caseServer_.getTime()
                    );
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
        }
    }

    //--------------------------------------------------------------------------

    private class RunCalculationAction
    extends AbstractAction
    {
        RunCalculationAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("RunCalculationImage")
            );
            putValue(Action.NAME, "Start Calculation");
            putValue(Action.SHORT_DESCRIPTION, "Start Calculation");
            putValue(Action.LONG_DESCRIPTION, "Start Calculation");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(App.getRootFrame());

            try
            {
                if (caseServer_ != null)
                {
                    boolean doRun =
                        checkAndSave
                        (
                            true,
                            caseServer_.caseName()
                          + " : case has been changed."
                          + " Save before running?",
                            "Save Case"
                        );
                    if (doRun)
                    {
                        RunPanel runPanel = new RunPanel
                        (
                            App.getRootFrame(),
                            caseBrowser_,
                            caseServer_.caseRoot(),
                            caseServer_.caseName(),
                            caseServer_.application().name(),
                            true,               // Is application
                            "controlDict",
                            null
                        );
                        runPanel.setVisible(true);
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

    // Not used anymore. Use runPanel always.
    private class StartCalculationAction
    extends AbstractAction
    {
        StartCalculationAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("StartCalculationImage")
            );
            putValue(Action.NAME, "Start Calculation Now");
            putValue(Action.SHORT_DESCRIPTION, "Start Calculation Now");
            putValue(Action.LONG_DESCRIPTION, "Start Calculation Now");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(App.getRootFrame());

            try
            {
                if (caseServer_ != null)
                {
                    boolean doRun = checkAndSave
                        (
                            true,
                            caseServer_.caseName() + " : case has been changed."
                            + " Save before running?",
                            "Save Case"
                        );
                    if (doRun)
                    {
                        // Start calculation.
                        int pid = caseServer_.runCase("");

                        System.out.println
                        ("Application started with pid " + pid);
                    }
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
        }
    }

    //--------------------------------------------------------------------------

    private class StopCalculationAction
    extends AbstractAction
    {
        StopCalculationAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("StopCalculationImage")
            );
            putValue(Action.NAME, "Stop Calculation");
            putValue(Action.SHORT_DESCRIPTION, "Stop Calculation");
            putValue(Action.LONG_DESCRIPTION, "Stop Calculation");
        }

        public void actionPerformed(ActionEvent e)
        {

            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(App.getRootFrame());

            try
            {
                if (caseServer_ != null)
                {
                    // Stop a calculation.
                    caseServer_.killCase();
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
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}





