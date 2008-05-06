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

import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.util.*;
import javax.swing.*;
import javax.swing.tree.*;
import javax.swing.event.*;
import java.awt.*;

import org.omg.CosNaming.*;
import org.omg.CORBA.StringHolder;

import FoamXServer.ApplicationDescriptor;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;
import FoamXServer.FoamXSYSError;
import FoamXServer.CaseBrowser.ICaseBrowser;
import FoamXServer.CaseDescriptor;
import FoamXServer.CaseServer.ICaseServer;
import FoamXServer.CaseServer.ICaseServerHolder;
import FoamXServer.CaseServer.IFoamProperties;
import FoamXServer.CasePostServer.ICasePostServer;
import FoamXServer.CasePostServer.ICasePostServerHolder;
import FoamXServer.IDictionaryEntry;
import FoamXServer.IDictionaryEntryHolder;
import FoamXServer.ApplicationDescriptor;

import FoamX.App;
import FoamX.Exceptions.FoamXException;
import FoamX.Util.BusyCursor;
import FoamX.Util.CollapseTreeAction;
import FoamX.Util.ExpandTreeAction;
import FoamX.Util.FoamXTreeRenderer;
import FoamX.Tools.InvokeUtilityAction;
import FoamX.ToolbarManagement.Toolbar;
import FoamX.ProcessManagement.ProcessEditor;
import FoamX.Editors.CompoundEditor;
import FoamX.Editors.DictionaryEntryEditor.EntryCache.DictionaryEntryCache;
import FoamX.Editors.DictionaryEntryEditor.EntryCache.DictionaryCache;
import FoamX.Editors.DictionaryEntryEditor.EntryCache.WordCache;
import FoamX.Editors.DictionaryEntryEditor.EntryCache.ListCache;

public class CaseBrowserPanel
    extends javax.swing.JPanel
    implements FoamX.Util.TabSelection
{
    //--------------------------------------------------------------------------

    public static int SELECT_HOST_MODE = 1;
    public static int SELECT_ROOT_MODE = 2;
    public static int SELECT_CASE_MODE = 3;
    public static int OPEN_CASE_MODE = 4;

    protected CaseBrowserModel model_;
    protected EventListenerList listenerList_;
    protected Toolbar toolBar_;
    protected int selectionMode_;

    // Case Browser actions.
    StartCaseBrowserAction startCaseBrowserAction_;
    StopCaseBrowserAction stopCaseBrowserAction_;
    StartProcessEditorAction startProcessEditorAction_;
    EditFoamControlDictAction editFoamControlDictAction_;
    OpenRootAction openRootAction_;
    OpenCaseAction openCaseAction_;
    //OpenCasePostAction openCasePostAction_;
    CreateCaseAction createCaseAction_;
    //ImportCaseAction importCaseAction_;
    DeleteCaseAction deleteCaseAction_;
    CloneCaseAction cloneCaseAction_;
    UnlockCaseAction unlockCaseAction_;

    SelectHostAction selectHostAction_;
    SelectRootAction selectRootAction_;
    SelectCaseAction selectCaseAction_;

    RefreshAction refreshAction_;
    ExpandTreeAction expandAllAction_;
    CollapseTreeAction collapseAllAction_;


    //--------------------------------------------------------------------------
    /** CaseBrowserPanel constructor.
      *
      * Implements a panel where users can select hosts to run caseBrowser on
      * and select root and case.
      * Works in one of three modes, determined by selectionMode argument:
      * - caseserver mode: normal operation. Double clicking on case starts
      *   editing of case. (i.e. fires caseStatus events)
      * - host selection mode: used for selecting hosts.
      *   Double clicking on hosts fires hostSelection events
      * - root selection mode: used for selecting root directories.
      *   Double clicking on roots fires rootSelection events
      * - case selection mode: used for selecting case. Double clicking on
      *   case just selects case. (i.e. fires caseSelection events)
      *
      * For selection mode one normally reuses the CaseBrowserModel
      * so only one ICaseBrowser per node runs.
      */
    public CaseBrowserPanel(CaseBrowserModel model, int selectionMode)
    {
        try
        {
            selectionMode_ = selectionMode;

            // Initialise the GUI components.
            initComponents();

            // Create listener list.
            listenerList_ = new EventListenerList();

            // Set the tree model and renderer.
            if (model != null)
            {
                model_ = model;
            }
            else
            {
                model_ = new CaseBrowserModel();
            }
            tree_.setModel(model_);
            tree_.setCellRenderer(new FoamXTreeRenderer());
            tree_.putClientProperty
            (
                "JTree.lineStyle",
                App.getResources().getResourceString("Tree.LineStyle")
            );


            // Initialise actions.
            startCaseBrowserAction_     = new StartCaseBrowserAction();
            stopCaseBrowserAction_      = new StopCaseBrowserAction();
            startProcessEditorAction_   = new StartProcessEditorAction();
            editFoamControlDictAction_  = new EditFoamControlDictAction();

            openRootAction_             = new OpenRootAction();
            openCaseAction_             = new OpenCaseAction();
            //openCasePostAction_         = new OpenCasePostAction();
            createCaseAction_           = new CreateCaseAction();
            //importCaseAction_           = new ImportCaseAction();
            deleteCaseAction_           = new DeleteCaseAction();
            cloneCaseAction_            = new CloneCaseAction();
            unlockCaseAction_           = new UnlockCaseAction();

            selectHostAction_           = new SelectHostAction();
            selectRootAction_           = new SelectRootAction();
            selectCaseAction_           = new SelectCaseAction();

            refreshAction_              = new RefreshAction();
            expandAllAction_            = new ExpandTreeAction(tree_);
            collapseAllAction_          = new CollapseTreeAction(tree_);

            if (selectionMode_ == OPEN_CASE_MODE)
            {
                // Initialise toolbar.
                toolBar_ = App.getToolbarManager().addToolbar
                (
                    "FoamXCaseBrowser",
                    "Case Browser"
                );

                // Initialise toolbar.
                addToolbarButton(startProcessEditorAction_);
                toolBar_.addSeparator();
                addToolbarButton(startCaseBrowserAction_);
                addToolbarButton(stopCaseBrowserAction_);
                toolBar_.addSeparator();
                addToolbarButton(createCaseAction_);
                //addToolbarButton(importCaseAction_);
                addToolbarButton(openCaseAction_);
                //addToolbarButton(openCasePostAction_);
                addToolbarButton(deleteCaseAction_);
                addToolbarButton(cloneCaseAction_);
                addToolbarButton(unlockCaseAction_);
                toolBar_.addSeparator();
                addToolbarButton(editFoamControlDictAction_);
                toolBar_.addSeparator();
                addToolbarButton(refreshAction_);
                addToolbarButton(expandAllAction_);
                addToolbarButton(collapseAllAction_);

            }
            // Reset all actions.
            resetActions(true);

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

            if (model == null)
            {
                // Make sure the tree is displayed properly.
                updateUI();
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Close all case browsers. */
    public void shutdown()
    {
        try
        {
            // Close all case browsers.
            model_.shutdown();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Get selection mode */
    public int getSelectionMode()
    {
        return selectionMode_;
    }

    //--------------------------------------------------------------------------
    /** Set selection mode */
    public void setSelectionMode(int selectionMode)
    {
        selectionMode_ = selectionMode;
    }

    //--------------------------------------------------------------------------
    /** Get tree model */
    public CaseBrowserModel getModel()
    {
        return model_;
    }

    //--------------------------------------------------------------------------
    /** Set tree model */
    public void setModel(CaseBrowserModel model)
    {
        model_ = model;
    }

    //--------------------------------------------------------------------------
    /** Refresh the case browser model. */
    public void refreshCaseBrowser(boolean forceReread)
    {
        //// Refresh the hosts lists
        //model_.refreshHosts();

        // Refresh the case browser tree model.
        model_.refresh(forceReread);
        tree_.updateUI();
    }

    //--------------------------------------------------------------------------
    /** Refresh the case browser model for single directory only */
    public void refreshRoot(String hostName, String caseRoot)
    {
        // Refresh the case browser tree model.
        model_.refreshRoot(hostName, caseRoot);
        tree_.updateUI();
    }

    //--------------------------------------------------------------------------
    /** Refresh the case browser model for single directory only */
    public void refreshRoot(ICaseBrowser caseBrowser, String caseRoot)
        throws FoamXError
    {
        // Get hostname from browser
        StringHolder holder = new StringHolder();
        caseBrowser.getHostName(holder);

        refreshRoot(holder.value, caseRoot);
    }

    //--------------------------------------------------------------------------
    /** Refresh a single case in the case browser model. */
    public void refreshCaseNode(String caseRoot, String caseName)
    {
        // Refresh the case browser tree model.
        model_.refreshCaseNode
        (
            caseRoot,
            caseName
        );
        tree_.updateUI();
    }

    //--------------------------------------------------------------------------
    /** Delete a single case in the case browser model. */
    public void deleteCaseNode(String caseRoot, String caseName)
    {
        // Refresh the case browser tree model.
        model_.deleteCaseNode
        (
            caseRoot,
            caseName
        );
        tree_.updateUI();
    }

    //--------------------------------------------------------------------------
    /** Refresh the case browser model for hosts as well as roots. */
    public void refreshAll()
    {
        // Refresh the hosts lists
        model_.refreshHosts();

        // Reread roots/cases
        model_.refresh(true);

        tree_.updateUI();
    }

    //--------------------------------------------------------------------------
    /** Set host browser status. */
    public void setHostOffLine(String hostName, String msg)
    {
        // Refresh the host browser with information on a certain host
        model_.setHostOffLine
        (
            hostName,
            msg
        );
        tree_.updateUI();
    }

    //--------------------------------------------------------------------------
    /** Handle commandline opening of various */
    public boolean serverStart
    (
        String hostName,
        String caseRoot,
        String caseName,
        String action
    )
    {
        if (hostName != null)
        {
            Rectangle rect = new Rectangle();
            int row = getHostPos(hostName, rect);
            if (row == -1)
            {
                return false;
            }

            // Insert host-open event
            OpenHostThread openHost = new OpenHostThread(this, hostName);
            SwingUtilities.invokeLater(openHost);
        }

        if (caseRoot != null)
        {
            // Insert root-open event
            OpenRootThread openRoot = new OpenRootThread(this, caseRoot);
            SwingUtilities.invokeLater(openRoot);
        }

        if ((caseName != null) && (action != null))
        {
            if (action.equals("open"))
            {
                // Insert case-open event
                OpenCaseThread openCase = new OpenCaseThread(this, caseName);
                SwingUtilities.invokeLater(openCase);
            }
            else if (action.equals("post"))
            {
                // Insert case-post event
                PostCaseThread postCase = new PostCaseThread(this, caseName);
                SwingUtilities.invokeLater(postCase);
            }
        }
        return true;
    }

    //--------------------------------------------------------------------------
    /** Helper class to convert tree entry into x,y */
    private int getHostPos(String hostName, Rectangle rect)
    {
        for(int row = 0; row < tree_.getRowCount(); row++)
        {
            Rectangle rowRect = tree_.getRowBounds(row);

            TreePath tp = tree_.getPathForLocation(rowRect.x, rowRect.y);
            if (tp != null)
            {
                // Get node info.
                DefaultMutableTreeNode hostNode = 
                    (DefaultMutableTreeNode)tp.getLastPathComponent();
                ContextInfo context = (ContextInfo)hostNode.getUserObject();
                String contextName = context.getText();

                if (context.isHost() && contextName.equals(hostName))
                {
                    // Make sure the node is selected.
                    tree_.getSelectionModel().clearSelection();
                    tree_.getSelectionModel().setSelectionPath(tp);

                    rect.x = rowRect.x;
                    rect.y = rowRect.y;

                    return row;
                }
            }
        }

        return -1;
    }

    //--------------------------------------------------------------------------
    /** Helper class to convert tree entry into x,y */
    private int getRootPos(String caseRoot, Rectangle rect)
    {
        for(int row = 0; row < tree_.getRowCount(); row++)
        {
            Rectangle rowRect = tree_.getRowBounds(row);

            TreePath tp = tree_.getPathForLocation(rowRect.x, rowRect.y);
            if (tp != null)
            {
                // Get node info.
                DefaultMutableTreeNode node = 
                    (DefaultMutableTreeNode)tp.getLastPathComponent();
                ContextInfo context = (ContextInfo)node.getUserObject();

                if
                (
                    context.isCaseRoot()
                 && (caseRoot.equals(context.getCaseRoot()))
                )
                {
                    // Make sure the node is selected.
                    tree_.getSelectionModel().clearSelection();
                    tree_.getSelectionModel().setSelectionPath(tp);

                    rect.x = rowRect.x;
                    rect.y = rowRect.y;

                    return row;
                }
            }
        }

        return -1;
    }

    //--------------------------------------------------------------------------
    /** Helper class to convert tree entry into x,y */
    private int getCaseNamePos(String caseName, Rectangle rect)
    {
        for(int row = 0; row < tree_.getRowCount(); row++)
        {
            Rectangle rowRect = tree_.getRowBounds(row);

            TreePath tp = tree_.getPathForLocation(rowRect.x, rowRect.y);
            if (tp != null)
            {
                // Get node info.
                DefaultMutableTreeNode node = 
                    (DefaultMutableTreeNode)tp.getLastPathComponent();
                ContextInfo context = (ContextInfo)node.getUserObject();

                if
                (
                    context.isCaseName()
                 && (context.getCaseName().equals(caseName))
                )
                {
                    // Make sure the node is selected.
                    tree_.getSelectionModel().clearSelection();
                    tree_.getSelectionModel().setSelectionPath(tp);

                    rect.x = rowRect.x;
                    rect.y = rowRect.y;

                    return row;
                }
            }
        }

        return -1;
    }

    //--------------------------------------------------------------------------
    /** Thread class for mouse selection action. */
    protected class OpenHostThread implements Runnable
    {
        JComponent comp_;
        String hostName_;

        public OpenHostThread(JComponent comp, String hostName)
        {
            super();
            comp_ = comp;
            hostName_ = hostName;
        }

        public void run()
        {
            Rectangle rect = new Rectangle();
            int row = getHostPos(hostName_, rect);

            if (row != -1)
            {
                // Found host at row. Do as if double-clicked on it.
                if (!tree_.isExpanded(row))
                {
                    MouseEvent evt =
                        new MouseEvent
                        (
                            comp_,          // component
                            0,              // id
                            8888888888888L, // when
                            0,              // modifiers
                            rect.x,         // x
                            rect.y,         // y
                            2,              // click count
                            false           // popup trigger
                        );
                    OnMouseClicked(evt);
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    /** Thread class for mouse selection action. */
    protected class OpenRootThread implements Runnable
    {
        JComponent comp_;
        String caseRoot_;

        public OpenRootThread(JComponent comp, String caseRoot)
        {
            super();
            comp_ = comp;
            caseRoot_ = caseRoot;
        }

        public void run()
        {
            Rectangle rect = new Rectangle();
            int row = getRootPos(caseRoot_, rect);

            if (row != -1)
            {
                // Found caseRoot at row. Expand tree to show cases.
                if (!tree_.isExpanded(row))
                {
                    openRootAction_.actionPerformed(null);
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    /** Thread class for mouse selection action. */
    protected class OpenCaseThread implements Runnable
    {
        JComponent comp_;
        String caseName_;

        public OpenCaseThread(JComponent comp, String caseName)
        {
            super();
            comp_ = comp;
            caseName_ = caseName;
        }

        public void run()
        {
            Rectangle rect = new Rectangle();
            int row = getCaseNamePos(caseName_, rect);

            if (row != -1)
            {
                openCaseAction_.actionPerformed(null);
            }
        }
    }

    //--------------------------------------------------------------------------
    /** Thread class for mouse selection action. */
    protected class PostCaseThread implements Runnable
    {
        JComponent comp_;
        String caseName_;

        public PostCaseThread(JComponent comp, String caseName)
        {
            super();
            comp_ = comp;
            caseName_ = caseName;
        }

        public void run()
        {
            Rectangle rect = new Rectangle();
            int row = getCaseNamePos(caseName_, rect);

            if (row != -1)
            {
                //openCasePostAction_.actionPerformed(null);
            }
        }
    }


    //--------------------------------------------------------------------------
    /** Reset all actions to the default (disabled) state. */
    protected void resetActions(boolean tabSelected)
    {
        startCaseBrowserAction_.setEnabled(false);
        stopCaseBrowserAction_.setEnabled(false);
        startProcessEditorAction_.setEnabled(false);
        editFoamControlDictAction_.setEnabled(false);
        openRootAction_.setEnabled(false);
        openCaseAction_.setEnabled(false);
        //openCasePostAction_.setEnabled(false);
        createCaseAction_.setEnabled(false);
        //importCaseAction_.setEnabled(false);
        deleteCaseAction_.setEnabled(false);
        cloneCaseAction_.setEnabled(false);
        unlockCaseAction_.setEnabled(false);

        selectHostAction_.setEnabled(false);
        selectRootAction_.setEnabled(false);
        selectCaseAction_.setEnabled(false);

        // Always available if this tab is selected.
        refreshAction_.setEnabled(tabSelected);
        expandAllAction_.setEnabled(tabSelected);
        collapseAllAction_.setEnabled(tabSelected);
    }

    //--------------------------------------------------------------------------
    /** Set all actions according to the current selection. */
    protected void setActions()
    {
        try
        {
            // Reset all actions.
            resetActions(true);

            TreePath tp = tree_.getSelectionModel().getSelectionPath();
            if (tp != null)
            {
                DefaultMutableTreeNode node =
                    (DefaultMutableTreeNode)tp.getLastPathComponent();
                ContextInfo context = (ContextInfo)node.getUserObject();

                if (context != null)
                {
                    if (context.isHost())
                    {
                        if (selectionMode_ == SELECT_HOST_MODE)
                        {
                            selectHostAction_.setEnabled(true);
                        }
                        else
                        {
                            if (context.getCaseBrowser() == null)
                            {
                                startCaseBrowserAction_.setEnabled(true);
                            }
                            else
                            {
                                stopCaseBrowserAction_.setEnabled(true);
                                startProcessEditorAction_.setEnabled(true);
                                editFoamControlDictAction_.setEnabled(true);
                            }
                        }
                    }
                    else if (context.isCaseRoot())
                    {
                        if (selectionMode_ == OPEN_CASE_MODE)
                        {
                            startProcessEditorAction_.setEnabled(true);
                            editFoamControlDictAction_.setEnabled(true);
                            createCaseAction_.setEnabled(true);
                            //importCaseAction_.setEnabled(true);
                            openRootAction_.setEnabled(true);
                        }
                        else if (selectionMode_ == SELECT_ROOT_MODE)
                        {
                            selectRootAction_.setEnabled(true);
                        }
                        else if (selectionMode_ == SELECT_CASE_MODE)
                        {
                            openRootAction_.setEnabled(true);
                        }
                    }
                    else if (context.isCaseName())
                    {
                        if (selectionMode_ == OPEN_CASE_MODE)
                        {
                            startProcessEditorAction_.setEnabled(true);
                            editFoamControlDictAction_.setEnabled(true);
                            openCaseAction_.setEnabled(true);
                            //openCasePostAction_.setEnabled(true);
                            deleteCaseAction_.setEnabled(true);
                            cloneCaseAction_.setEnabled(true);
                            //if (!context.isCaseManaged())
                            //{
                            //    importCaseAction_.setEnabled(true);
                            //}
                            unlockCaseAction_.setEnabled(true);
                        }
                        else if (selectionMode_ == SELECT_CASE_MODE)
                        {
                            selectCaseAction_.setEnabled(true);
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

    //--------------------------------------------------------------------------

    protected void addToolbarButton(AbstractAction action)
    {
        // Add new button to the toolbar using the action object.
        javax.swing.JButton button = toolBar_.add(action);

        // Set the tooltip text.
        button.setToolTipText
        (
            (String)action.getValue(Action.SHORT_DESCRIPTION)
        );
        button.setText("");
        button.setFont(toolBar_.getFont());      // Use same font as toolbar.
    }

    //--------------------------------------------------------------------------

    protected void buildUtilsMenu(JPopupMenu contextMenu, ContextInfo context)
    {
        try
        {
            // Get case browser reference from context info object.
            ICaseBrowser caseBrowser = context.getCaseBrowser();
            if (caseBrowser != null)
            {
                boolean caseNode = context.isCaseName();

                Hashtable utilityMenuMap = new Hashtable(5);

                // Construct root utilities menu and add to context menu.
                JMenu utilitiesMenu = new JMenu("Foam Utilities");
                utilitiesMenu.setFont(contextMenu_.getFont());
                contextMenu_.add(utilitiesMenu);

                // Get utility descriptors.
                ApplicationDescriptor[] utilities =
                    caseBrowser.foamProperties().utilities();

                for (int i = 0; i <utilities.length; i++)
                {
                    // Only add case utilities to case nodes.
                    // Add non case utility anyway.
//                    if
//                    (
//                        !utilities[i].caseTool
//                     || (utilities[i].caseTool && caseNode)
//                    )
//                    {
                    // Build category sub-menus.
                    JMenu menu = utilitiesMenu;

                    if
                    (
                        utilities[i].category != null
                     && utilities[i].category.length() > 0
                    )
                    {
                        // Add all required sub-menus.
                        StringTokenizer tok= new StringTokenizer
                        (
                            utilities[i].category,
                            "/"
                        );
                        String catKey = new String();
                        while (tok.hasMoreTokens())
                        {
                            String category = tok.nextToken();
                            catKey = catKey + category;
                            if (utilityMenuMap.containsKey(catKey))
                            {
                                menu = (JMenu)utilityMenuMap.get(catKey);
                            }
                            else
                            {
                                JMenu catMenu = new JMenu(category);
                                catMenu.setFont(contextMenu_.getFont());
                                menu.add(catMenu);
                                utilityMenuMap.put(catKey, catMenu);
                                menu = catMenu;
                            }
                        }
                    }

                    // Add a menu item for this utility.
                    JMenuItem utilityItem = new JMenuItem(utilities[i].name);
                    utilityItem.setFont(contextMenu_.getFont());

                    // Hook up to invoke utility action.
                    utilityItem.addActionListener
                    (
                        new InvokeUtilityAction
                        (
                            caseBrowser,
                            context.getCaseRoot(),
                            context.getCaseName()
                        )
                    );

                    // Specify which utility to invoke and on which case
                    // browser.
                    utilityItem.putClientProperty
                    (
                        "utilityDescriptor",
                        utilities[i]
                    );

                    // Add to menu.
                    menu.add(utilityItem);
//                    }
                }
            }
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
        setActions();
    }

    public void tabDeselected()
    {
        resetActions(false);
    }

    //--------------------------------------------------------------------------
    //---- CaseStatusListener Event Source Methods
    //--------------------------------------------------------------------------

    public synchronized void addCaseStatusListener(CaseStatusListener l)
    {
        listenerList_.add(CaseStatusListener.class, l);
    }

    //--------------------------------------------------------------------------

    public synchronized void removeCaseStatusListener(CaseStatusListener l)
    {
        listenerList_.remove(CaseStatusListener.class, l);
    }

    //--------------------------------------------------------------------------

    protected synchronized void fireCaseOpened
    (
        String caseRoot,
        String caseName,
        ICaseServer caseServer,
        ICaseBrowser caseBrowser
    )
    {
        // Create event object.
        CaseStatusEvent evt = new CaseStatusEvent
        (
            this,
            caseRoot,
            caseName,
            caseServer,
            null,
            caseBrowser
        );

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseStatusListener.class)
            {
                ((CaseStatusListener)listeners[i+1]).caseOpened(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected synchronized void fireCasePostOpened
    (
        String caseRoot,
        String caseName,
        ICasePostServer casePostServer,
        ICaseBrowser caseBrowser
    )
    {
        // Create event object.
        CaseStatusEvent evt = new CaseStatusEvent
        (
            this,
            caseRoot,
            caseName,
            null,
            casePostServer,
            caseBrowser
        );

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseStatusListener.class)
            {
                ((CaseStatusListener)listeners[i+1]).casePostOpened(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected synchronized void fireCaseClosed(String caseRoot, String caseName)
    {
        // Create event object.
        CaseStatusEvent evt = new CaseStatusEvent(this, caseRoot, caseName);

        // Process the listeners last to first, notifying those that are
        // interested in this event.
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

    protected synchronized void fireCaseDeleted
    (
        String caseRoot,
        String caseName
    )
    {
        // Create event object.
        CaseStatusEvent evt = new CaseStatusEvent(this, caseRoot, caseName);

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseStatusListener.class)
            {
                ((CaseStatusListener)listeners[i+1]).caseDeleted(evt);
            }
        }
    }


    //--------------------------------------------------------------------------
    //---- CaseSelectionListener Event Source Methods
    //--------------------------------------------------------------------------

    public synchronized void addCaseSelectionListener(CaseSelectionListener l)
    {
        listenerList_.add(CaseSelectionListener.class, l);
    }

    //--------------------------------------------------------------------------

    public synchronized void removeCaseSelectionListener(CaseSelectionListener l)
    {
        listenerList_.remove(CaseSelectionListener.class, l);
    }

    //--------------------------------------------------------------------------

    protected synchronized void fireCaseSelected
    (
        String hostName,
        String caseRoot,
        String caseName,
        ICaseBrowser caseBrowser
    )
    {
        // Create event object.
        CaseSelectionEvent evt =
            new CaseSelectionEvent
            (
                this,
                hostName, 
                caseRoot, 
                caseName,
                caseBrowser
            );

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseSelectionListener.class)
            {
                ((CaseSelectionListener)listeners[i+1]).caseSelected(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected synchronized void fireRootSelected
    (
        String hostName,
        String caseRoot,
        ICaseBrowser caseBrowser
    )
    {
        // Create event object.
        CaseSelectionEvent evt =
            new CaseSelectionEvent
            (
                this,
                hostName, 
                caseRoot, 
                null,
                caseBrowser
            );

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseSelectionListener.class)
            {
                ((CaseSelectionListener)listeners[i+1]).rootSelected(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected synchronized void fireHostSelected(String hostName)
    {
        // Create event object.
        CaseSelectionEvent evt =
            new CaseSelectionEvent
            (
                this,
                hostName, 
                null, 
                null,
                null
            );

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseSelectionListener.class)
            {
                ((CaseSelectionListener)listeners[i+1]).hostSelected(evt);
            }
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
        contextMenu_ = new javax.swing.JPopupMenu();
        scrollPane_ = new javax.swing.JScrollPane();
        tree_ = new javax.swing.JTree();
        
        contextMenu_.setFont(new java.awt.Font("Dialog", 0, 10));
        
        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        tree_.setFont(new java.awt.Font("Dialog", 0, 10));
        tree_.setShowsRootHandles(true);
        tree_.addMouseListener(new java.awt.event.MouseAdapter()
        {
            public void mouseClicked(java.awt.event.MouseEvent evt)
            {
                OnMouseClicked(evt);
            }
        });
        
        tree_.addTreeExpansionListener(new javax.swing.event.TreeExpansionListener()
        {
            public void treeExpanded(javax.swing.event.TreeExpansionEvent evt)
            {
                OnTreeExpand(evt);
            }
            public void treeCollapsed(javax.swing.event.TreeExpansionEvent evt)
            {
            }
        });
        
        tree_.addTreeWillExpandListener(new javax.swing.event.TreeWillExpandListener()
        {
            public void treeWillExpand(javax.swing.event.TreeExpansionEvent evt)
            throws javax.swing.tree.ExpandVetoException
            {
                OnTreeWillExpand(evt);
            }
            public void treeWillCollapse(javax.swing.event.TreeExpansionEvent evt)
            throws javax.swing.tree.ExpandVetoException
            {
            }
        });
        
        scrollPane_.setViewportView(tree_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.SOUTHWEST;
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        add(scrollPane_, gridBagConstraints1);
        
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

  private void OnTreeWillExpand(javax.swing.event.TreeExpansionEvent evt)
    {//GEN-FIRST:event_OnTreeWillExpand
        // Add your handling code here:
  }//GEN-LAST:event_OnTreeWillExpand

    //--------------------------------------------------------------------------

  private void OnTreeExpand(javax.swing.event.TreeExpansionEvent evt)
    {//GEN-FIRST:event_OnTreeExpand
        // Add your handling code here:
    // Get the node object for the current selection.
//    TreePath tp = evt.getPath();
//    if (tp != null)
//    {
//        // Make sure the node is selected.
//        tree_.getSelectionModel().clearSelection();
//        tree_.getSelectionModel().setSelectionPath(tp);
//
//        DefaultMutableTreeNode node =
//            (DefaultMutableTreeNode)tp.getLastPathComponent();
//        ContextInfo context = (ContextInfo)node.getUserObject();
//
//        if (context.isCaseRoot())
//        {
//
//    System.out.println("OnTreeExpand root:" + context.getCaseRoot());
//
//            DefaultMutableTreeNode parentNode =
//                (DefaultMutableTreeNode)
//                tp.getParentPath().getLastPathComponent();
//            ContextInfo parentContext =
//                (ContextInfo)parentNode.getUserObject();
//
//    System.out.println("OnTreeExpand host:" + parentContext.getText());
//
//            if (parentContext.isHost())
//            {
//                String rootDir = context.getCaseRoot();
//                String hostName = parentContext.getText();
//                refreshRoot(hostName, rootDir);
//
//                tree_.expandPath(tp);
//            }
//        }
//    }
  }//GEN-LAST:event_OnTreeExpand

    //--------------------------------------------------------------------------

  private void OnMouseClicked (java.awt.event.MouseEvent evt)
    {//GEN-FIRST:event_OnMouseClicked
     if (selectionMode_ == OPEN_CASE_MODE)
     {
         handleMouseOpen(evt);
     }
     else if 
     (
         (selectionMode_ == SELECT_HOST_MODE)
      || (selectionMode_ == SELECT_ROOT_MODE)
      || (selectionMode_ == SELECT_CASE_MODE)
     )
     {
         handleMouseSelection(evt);
     }
     else
     {
         System.out.println("ILLEGAL:CaseBrowserPanel::OnMouseClicked");
     }
  }//GEN-LAST:event_OnMouseClicked
  

    //--------------------------------------------------------------------------

   
    private void handleMouseOpen(java.awt.event.MouseEvent evt)
    {
        try
        {
            // Get the node object for the current selection.
            TreePath tp = tree_.getPathForLocation(evt.getX(), evt.getY());
            if (tp != null)
            {
                // Make sure the node is selected.
                tree_.getSelectionModel().clearSelection();
                tree_.getSelectionModel().setSelectionPath(tp);

                DefaultMutableTreeNode node =
                    (DefaultMutableTreeNode)tp.getLastPathComponent();
                ContextInfo context = (ContextInfo)node.getUserObject();

                // Check for popup menu trigger.
                //if (evt.isPopupTrigger())
                // Check explicitly for right mouse button.
                if (evt.getModifiers() == evt.BUTTON3_MASK)
                {
                    // Construct an appropriate context menu,
                    contextMenu_.removeAll();
                    java.awt.Font font = contextMenu_.getFont();

                    if (context == null)
                    {
                        // Do not have a licence manager reference.
                        // Can only refresh the tree.
                        contextMenu_.add(refreshAction_).setFont(font);
                    }
                    else if (context.isRootContext())
                    {
                        // Can only refresh or expand/collapse the tree.
                        contextMenu_.add(refreshAction_).setFont(font);
                        contextMenu_.add(expandAllAction_).setFont(font);
                        contextMenu_.add(collapseAllAction_).setFont(font);
                    }
                    else if (context.isHost())
                    {
                        // Host node
                        if (context.getCaseBrowser() == null)
                        {
                            // Allow start caseBrowser if non running
                            contextMenu_.add(startCaseBrowserAction_).setFont
                            (
                                font
                            );
                        }
                        else
                        {
                            // Allow various and utilities invocation
                            contextMenu_.add(stopCaseBrowserAction_).setFont
                            (
                                font
                            );
                            contextMenu_.add(createCaseAction_).setFont(font);
                            //contextMenu_.add(importCaseAction_).setFont(font);
                            contextMenu_.add(startProcessEditorAction_).setFont
                            (
                                font
                            );
                            contextMenu_.add(editFoamControlDictAction_).setFont
                            (
                                font
                            );
                            // Allow invoke utilities
                            buildUtilsMenu(contextMenu_, context);
                        }
                    }
                    else if (context.isCaseRoot())
                    {
                        // caseroot: open directory, create/import cases
                        contextMenu_.add(openRootAction_).setFont(font);
                        contextMenu_.add(createCaseAction_).setFont(font);
                        //contextMenu_.add(importCaseAction_).setFont(font);

                        // Allow processEditing on rootDir
                        contextMenu_.add(startProcessEditorAction_).setFont
                        (
                            font
                        );
                        contextMenu_.add(editFoamControlDictAction_).setFont
                        (
                            font
                        );
                        // Allow invoke utilities
                        buildUtilsMenu(contextMenu_, context);
                    }
                    else if (context.isCaseName())
                    {
                        // Can open, delete or unlock specific cases.
                        contextMenu_.add(openCaseAction_).setFont(font);
                        //contextMenu_.add(openCasePostAction_).setFont(font);
                        contextMenu_.add(deleteCaseAction_).setFont(font);
                        contextMenu_.add(cloneCaseAction_).setFont(font);
                        //contextMenu_.add(importCaseAction_).setFont(font);
                        contextMenu_.add(unlockCaseAction_).setFont(font);

                        // caseName: allow processEditing
                        contextMenu_.add(startProcessEditorAction_).setFont
                        (
                            font
                        );
                        // Allow invoke utilities
                        buildUtilsMenu(contextMenu_, context);
                    }

                    // Show context menu.
                    contextMenu_.show(tree_, evt.getX(), evt.getY());
                }
                else if (evt.getClickCount() >= 2)
                {
                    // Start a case browser if double clicking on a host node.
                    if (context.isHost() && context.getCaseBrowser() == null)
                    {
                        startCaseBrowserAction_.actionPerformed(null);
                    }
                    // Invoke the open case action if the user
                    // has double clicked on a case.
                    else if (context.isCaseName())
                    {
                        openCaseAction_.actionPerformed(null);
                    }
                    // Invoke the open root action if the user
                    // has double clicked on a root.
                    else if (context.isCaseRoot())
                    {
                        openRootAction_.actionPerformed(null);
                    }
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }


    //--------------------------------------------------------------------------

   
    //--------------------------------------------------------------------------
    /** Handle mouse clicks if in selection mode. */
    private void handleMouseSelection(java.awt.event.MouseEvent evt)
    {
        try
        {
            // Get the node object for the current selection.
            TreePath tp = tree_.getPathForLocation(evt.getX(), evt.getY());
            if (tp != null)
            {
                // Make sure the node is selected.
                tree_.getSelectionModel().clearSelection();
                tree_.getSelectionModel().setSelectionPath(tp);

                DefaultMutableTreeNode node =
                    (DefaultMutableTreeNode)tp.getLastPathComponent();
                ContextInfo context = (ContextInfo)node.getUserObject();

                // Check explicitly for right mouse button.
                if (evt.getModifiers() == evt.BUTTON3_MASK)
                {
                    // Construct an appropriate context menu,
                    contextMenu_.removeAll();
                    java.awt.Font font = contextMenu_.getFont();

                    if (context == null)
                    {
                        // Do not have a licence manager reference.
                        // Can only refresh the tree.
                        contextMenu_.add(refreshAction_).setFont(font);
                    }
                    else if (context.isRootContext())
                    {
                        // Can only refresh or expand/collapse the tree.
                        contextMenu_.add(refreshAction_).setFont(font);
                        contextMenu_.add(expandAllAction_).setFont(font);
                        contextMenu_.add(collapseAllAction_).setFont(font);
                    }
                    else if (context.isHost())
                    {
                        if (selectionMode_ == SELECT_HOST_MODE)
                        {
                            contextMenu_.add(selectHostAction_).setFont(font);
                        }
                        else if (context.getCaseBrowser() == null)
                        {
                            // Can start a case browser.
                            contextMenu_.add(startCaseBrowserAction_).setFont
                            (
                                font
                            );
                        }
                        else
                        {
                            // Can stop a case browser.
                            contextMenu_.add(stopCaseBrowserAction_).setFont
                            (
                                font
                            );
                        }
                    }
                    else if (context.isCaseRoot())
                    {
                        if (selectionMode_ == SELECT_ROOT_MODE)
                        {
                            contextMenu_.add(selectRootAction_).setFont(font);
                        }
                        else if (selectionMode_ == SELECT_CASE_MODE)
                        {
                            contextMenu_.add(openRootAction_).setFont(font);
                        }
                    }
                    else if (context.isCaseName())
                    {
                        if (selectionMode_ == SELECT_CASE_MODE)
                        {
                            contextMenu_.add(selectCaseAction_).setFont(font);
                        }
                    }

                    // Show context menu.
                    contextMenu_.show(tree_, evt.getX(), evt.getY());
                }
                else if (evt.getClickCount()>= 2)
                {
                    if (selectionMode_ == SELECT_HOST_MODE)
                    {
                        selectHostAction_.actionPerformed(null);
                    }
                    else if
                    (
                        context.isHost() && context.getCaseBrowser() == null
                    )
                    {
                        // Start a case browser if double clicking on a host
                        // node.
                        startCaseBrowserAction_.actionPerformed(null);
                    }
                    // Invoke the root/case selection if the user
                    // has double clicked on a root/case.
                    else if (context.isCaseRoot())
                    {
                        if (selectionMode_ == SELECT_ROOT_MODE)
                        {
                            selectRootAction_.actionPerformed(null);
                        }
                        else if (selectionMode_ == SELECT_CASE_MODE)
                        {
                            openRootAction_.actionPerformed(null);
                        }
                    }
                    else if (context.isCaseName())
                    {
                        if (selectionMode_ == SELECT_CASE_MODE)
                        {
                            selectCaseAction_.actionPerformed(null);
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

    //--------------------------------------------------------------------------

    protected void OnTreeSelectionChanged(TreeSelectionEvent evt)
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
    private javax.swing.JScrollPane scrollPane_;
    private javax.swing.JTree tree_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //---- Action Classes. CaseBrowser specific
    //--------------------------------------------------------------------------

    private class StartCaseBrowserAction
    extends AbstractAction
    {
        StartCaseBrowserAction()
        {
            putValue
            (
                Action.SMALL_ICON, App.getResources().getIcon
                (
                    "CaseBrowser.StartBrowserImage"
                )
            );
            putValue(Action.NAME, "Open Case Browser");
            putValue(Action.SHORT_DESCRIPTION, "Open Case Browser");
            putValue(Action.LONG_DESCRIPTION, "Open Case Browser");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            //BusyCursor busyCursor = new BusyCursor(CaseBrowserPanel.this);
            BusyCursor cursor = new BusyCursor(App.getRootFrame());

            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode hostNode = 
                        (DefaultMutableTreeNode)tp.getLastPathComponent();
                    ContextInfo context = (ContextInfo)hostNode.getUserObject();

                    // Check node type.
                    if
                    (
                        context.isHost()
                     && context.getCaseBrowser() == null
                     && model_.startCaseBrowser(context.getText())
                    )
                    {
                        // Expand the host node.
                        tree_.expandPath(tp);
                        tree_.updateUI();
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

    private class StopCaseBrowserAction
    extends AbstractAction
    {
        StopCaseBrowserAction()
        {
            putValue
            (
                Action.SMALL_ICON, App.getResources().getIcon
                (
                    "CaseBrowser.StopBrowserImage"
                )
            );
            putValue(Action.NAME, "Close Case Browser");
            putValue(Action.SHORT_DESCRIPTION, "Close Case Browser");
            putValue(Action.LONG_DESCRIPTION, "Close Case Browser");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(CaseBrowserPanel.this);

            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode hostNode = 
                        (DefaultMutableTreeNode)tp.getLastPathComponent();
                    ContextInfo context = (ContextInfo)hostNode.getUserObject();

                    // Check node type.
                    if (context.isHost() && context.getCaseBrowser() != null)
                    {
                        model_.stopCaseBrowser(context.getText());
                        tree_.updateUI();
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

    private class StartProcessEditorAction
    extends AbstractAction
    {
        StartProcessEditorAction()
        {
            putValue
            (
                Action.SMALL_ICON, App.getResources().getIcon
                (
                    "CaseBrowser.StartProcessEditorImage"
                )
            );
            putValue(Action.NAME, "Start Process Editor");
            putValue(Action.SHORT_DESCRIPTION, "Start Process Editor");
            putValue(Action.LONG_DESCRIPTION, "Start Process Editor");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(CaseBrowserPanel.this);

            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode hostNode = 
                        (DefaultMutableTreeNode)tp.getLastPathComponent();
                    ContextInfo context = (ContextInfo)hostNode.getUserObject();

                    // Check node type.
                    if (context.getCaseBrowser() != null)
                    {
                        // Start process editor as free standing window.
                        ProcessEditor processEditor =
                            new ProcessEditor
                            (
                                null,
                                context.getCaseBrowser()
                            );
                        processEditor.setVisible(true);
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

    private class EditFoamControlDictAction
    extends AbstractAction
    {
        EditFoamControlDictAction()
        {
            putValue
            (
                Action.SMALL_ICON, App.getResources().getIcon
                (
                    "CaseBrowser.EditFoamControlDictImage"
                )
            );
            putValue(Action.NAME, "Edit OpenFOAM controlDict");
            putValue(Action.SHORT_DESCRIPTION, "Edit OpenFOAM controlDict");
            putValue(Action.LONG_DESCRIPTION, "Edit OpenFOAM controlDict");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(CaseBrowserPanel.this);

            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode hostNode = 
                        (DefaultMutableTreeNode)tp.getLastPathComponent();
                    ContextInfo context = (ContextInfo)hostNode.getUserObject();

                    // Check node type.
                    if (context.getCaseBrowser() != null)
                    {
                        // We have a case browser. Load the system controlDict
                        // and start editor on it.
                        ICaseBrowser caseBrowser = context.getCaseBrowser();

                        IDictionaryEntryHolder dictHolder =
                            new IDictionaryEntryHolder();
                        caseBrowser.foamProperties().getFoamControlDict
                        (
                            dictHolder
                        );

                        if (dictHolder.value == null)
                        {
                            throw new FoamXException
                            (
                                "Could not find OpenFOAM controlDict"
                            );
                        }

                        CompoundEditor compoundEditor =
                            new CompoundEditor
                            (
                                App.getRootFrame(),
                                dictHolder.value
                            );
                        compoundEditor.setTitle("OpenFOAM Dictionary");
                        compoundEditor.setVisible(true);

                        if (dictHolder.value.modified())
                        {
                            System.out.println
                            (
                                "Saving modified OpenFOAM Dictionary"
                            );
                            dictHolder.value.save();
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

    // Connects to existing/starts new caseServer. Returns true if started.
    static private void getCaseServerReference
    (
        String caseRoot,
        String caseName,
        CaseDescriptor caseDescriptor,
        ICaseBrowser caseBrowser,
        ICaseServerHolder caseHolder
    ) throws FoamXException, FoamXSYSError, FoamXError, FoamXIOError
    {
        boolean caseOK = false;

        // First try to connect to already running server
        caseOK = caseBrowser.getCaseServerReference
        (
            caseRoot,
            caseName,
            caseHolder
        );
        if (caseOK)
        {
            //ICaseServer caseServer = caseHolder.value;
            String[] options = { "Use", "Cancel" };
            int ret =
                JOptionPane.showOptionDialog
                (
                    App.getRootFrame(),
                    "Found running Case Server"
                  + "\nUse it ?",
                    "Case Server",
                    JOptionPane.DEFAULT_OPTION,
                    JOptionPane.QUESTION_MESSAGE,
                    null,
                    options,
                    options[0]
                );
            if (ret == 0) return;
        }

        // Start case server
        if (caseDescriptor != null)
        {
            caseBrowser.openCase(caseDescriptor);
        }

        // Go into loop to detect when it has registered.
        int sleep = Integer.parseInt
        (
            App.getOptions().getProperty("FoamX.Sleep", "500")
        );
        int nRetries = Integer.parseInt
        (
            App.getOptions().getProperty("FoamX.NRetries", "20")
        );

        for (;;)
        {
            for (int i=0; i<nRetries; i++)
            {
                try
                {
                    Thread.currentThread().sleep(sleep);
                }
                catch (InterruptedException ie)
                {}

                caseOK = caseBrowser.getCaseServerReference
                (
                    caseRoot,
                    caseName,
                    caseHolder
                );

                if (caseOK) return;
            }

            if (!caseOK)
            {
                String[] options = { "Retry", "Exit" };
                int ret =
                    JOptionPane.showOptionDialog
                    (
                        App.getRootFrame(),
                        "Case Server cannot be contacted.",
                        "Case Server",
                        JOptionPane.DEFAULT_OPTION,
                        JOptionPane.ERROR_MESSAGE,
                        null,
                        options,
                        options[0]
                    );
                if (ret != 0)
                {
                    throw new FoamXException
                    (
                        "No caseServer started for "
                      + caseRoot + "/" + caseName
                    );
                }
            }
        }
    }

    //--------------------------------------------------------------------------

    static private void getCasePostServerReference
    (
        String caseRoot,
        String caseName,
        CaseDescriptor caseDescriptor,
        int nProcs,
        ICaseBrowser caseBrowser,
        ICasePostServerHolder caseHolder
    ) throws FoamXException, FoamXSYSError, FoamXError, FoamXIOError
    {
        boolean caseOK = false;

        // First try to connect to already running server
        caseOK = caseBrowser.getCasePostServerReference
        (
            caseRoot,
            caseName,
            nProcs,
            caseHolder
        );
        if (caseOK)
        {
            ICasePostServer casePostServer = caseHolder.value;
            if (casePostServer.nProcs() == nProcs)
            {
                String[] options = { "Use", "Cancel" };
                int ret =
                    JOptionPane.showOptionDialog
                    (
                        App.getRootFrame(),
                        "Found running Case Post Server with "
                      + casePostServer.nProcs() + " processors."
                      + "\nUse it ?",
                        "Case Post Server",
                        JOptionPane.DEFAULT_OPTION,
                        JOptionPane.QUESTION_MESSAGE,
                        null,
                        options,
                        options[0]
                    );
                if (ret == 0) return;
            }
            else
            {
                // Running but not for correct number of processors
                String[] options = { "Start", "Cancel" };
                int ret =
                    JOptionPane.showOptionDialog
                    (
                        App.getRootFrame(),
                        "Case Post Server running but with "
                      + casePostServer.nProcs() + " processors."
                      + "\nStart new casePostServer ?",
                        "Case Post Server",
                        JOptionPane.DEFAULT_OPTION,
                        JOptionPane.QUESTION_MESSAGE,
                        null,
                        options,
                        options[0]
                    );
                if (ret != 0)
                {
                    throw new FoamXException
                    (
                        "No casePostServer started for "
                      + caseRoot + "/" + caseName
                    );
                }
            }
        }

        // Start post server
        caseBrowser.openCasePost(caseDescriptor, nProcs);

        // Go into loop to detect when it has registered.
        int sleep = Integer.parseInt
        (
            App.getOptions().getProperty("FoamX.Sleep", "500")
        );
        int nRetries = Integer.parseInt
        (
            App.getOptions().getProperty("FoamX.NRetries", "20")
        );
        for (;;)
        {
            for (int i=0; i<nRetries; i++)
            {
                try
                {
                    Thread.currentThread().sleep(sleep);
                }
                catch (InterruptedException ie)
                {}

                caseOK = caseBrowser.getCasePostServerReference
                (
                    caseRoot,
                    caseName,
                    nProcs,
                    caseHolder
                );

                if (caseOK) return;
            }

            if (!caseOK)
            {
                String[] options = { "Retry", "Exit" };
                int ret =
                    JOptionPane.showOptionDialog
                    (
                        App.getRootFrame(),
                        "Case Post Server cannot be contacted.",
                        "Case Post Server",
                        JOptionPane.DEFAULT_OPTION,
                        JOptionPane.ERROR_MESSAGE,
                        null,
                        options,
                        options[0]
                    );
                if (ret != 0)
                {
                    throw new FoamXException
                    (
                        "No casePostServer started for "
                      + caseRoot + "/" + caseName
                    );
                }
            }
        }
    }

    //--------------------------------------------------------------------------

    static boolean openParallelDialog(int nProcs)
    {
        String[] options = { "Non-parallel", "Parallel" };
        int ret =
            JOptionPane.showOptionDialog
            (
                App.getRootFrame(),
                "Case is decomposed (" + nProcs + " processors)\n"
              + "Postprocess non-parallel or parallel.",
                "Case Post Server",
                JOptionPane.DEFAULT_OPTION,
                JOptionPane.ERROR_MESSAGE,
                null,
                options,
                options[1]
            );
        if (ret == 0)
        {
            return false;
        }
        else
        {
            return true;
        }
    }


    //--------------------------------------------------------------------------
    //---- Action Classes. Case specific
    //--------------------------------------------------------------------------

    private class OpenRootAction
        extends AbstractAction
    {
        OpenRootAction()
        {
            putValue
            (
                Action.SMALL_ICON, App.getResources().getIcon
                (
                    "CaseBrowser.OpenCaseImage"
                )
            );
            putValue(Action.NAME, "Open Root");
            putValue(Action.SHORT_DESCRIPTION, "Open Root");
            putValue(Action.LONG_DESCRIPTION, "Open Root");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(CaseBrowserPanel.this);

            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();

                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode node =
                        (DefaultMutableTreeNode)tp.getLastPathComponent(); 

                    ContextInfo context = 
                        (ContextInfo)node.getUserObject();

                    if (context.isCaseRoot())
                    {
                        DefaultMutableTreeNode parentNode =
                            (DefaultMutableTreeNode)
                            tp.getParentPath().getLastPathComponent();
                        ContextInfo parentContext =
                            (ContextInfo)parentNode.getUserObject();

                        if (parentContext.isHost())
                        {
                            String rootDir = context.getCaseRoot();
                            String hostName = parentContext.getText();
                            refreshRoot(hostName, rootDir);

                            tree_.expandPath(tp);
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

    private class OpenCaseAction
        extends AbstractAction
    {
        OpenCaseAction()
        {
            putValue
            (
                Action.SMALL_ICON, App.getResources().getIcon
                (
                    "CaseBrowser.OpenCaseImage"
                )
            );
            putValue(Action.NAME, "Open Case");
            putValue(Action.SHORT_DESCRIPTION, "Open Case");
            putValue(Action.LONG_DESCRIPTION, "Open Case");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            //MJ.No busy cursor. Upsets running from command line.
            // BusyCursor busyCursor = new BusyCursor(CaseBrowserPanel.this);

            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode nodeItem =
                        (DefaultMutableTreeNode)tp.getLastPathComponent(); 

                    ContextInfo contextInfo = 
                        (ContextInfo)nodeItem.getUserObject();

                    // Check node type.
                    if (contextInfo.isCaseName())
                    {
                        // Get case browser reference from context info object.
                        ICaseBrowser caseBrowser = contextInfo.getCaseBrowser();
                        CaseDescriptor caseDescriptor =
                            contextInfo.getCaseDescriptor();

                        // Final check to see if the case is locked.
                        if (caseBrowser.caseLocked(caseDescriptor))
                        {
                            // Issue warning.
                            JOptionPane.showMessageDialog
                            (
                                App.getRootFrame(),
                                "This case is locked."
                              + "\n(it is probably already being preprocessed)",
                                "Case Locked Error",
                                JOptionPane.ERROR_MESSAGE
                            );

                            // Mark as locked.
                            contextInfo.setCaseLocked(true);

                            // Exit stage left.
                            return;
                        }

                        // Open the selected case.
                        ICaseServerHolder caseHolder = new ICaseServerHolder();
                        try
                        {
                            getCaseServerReference
                            (
                                contextInfo.getCaseRoot(),
                                contextInfo.getCaseName(),
                                caseDescriptor,
                                caseBrowser,
                                caseHolder
                            );

                            // Update status from C++ side.
                            refreshCaseNode(contextInfo.getCaseRoot(), contextInfo.getCaseName());
                        }
                        catch (FoamXIOError ioErr)
                        {
                            // Mark the case as in error.
                            caseBrowser.caseIsInError(caseDescriptor);
                            refreshCaseNode(contextInfo.getCaseRoot(), contextInfo.getCaseName());

                            ioErr.errorMessage =
                                "The Open Case operation failed due to\n"
                              + ioErr.errorMessage;
                            throw ioErr;
                        }
                        catch (FoamXError fxErr)
                        {
                            // Mark the case as in error.
                            caseBrowser.caseIsInError(caseDescriptor);
                            refreshCaseNode(contextInfo.getCaseRoot(), contextInfo.getCaseName());

                            fxErr.errorMessage =
                                "The Open Case operation failed due to\n"
                              + fxErr.errorMessage;
                            throw fxErr;
                        }

                        ICaseServer caseServer = caseHolder.value;

                        // Fire case opened event.
                        fireCaseOpened
                        (
                            contextInfo.getCaseRoot(),
                            contextInfo.getCaseName(),
                            caseServer,
                            caseBrowser
                        );

                        if (!caseServer.managed())
                        {
                            // Mark the case as in error.
                            // (both locally and in the caseBrowser)
                            caseBrowser.caseIsInError(caseDescriptor);
                            refreshCaseNode(contextInfo.getCaseRoot(), contextInfo.getCaseName());

                            App.printMessage
                            (
                                App.DEBUGLEVEL_VERBOSE,
                                "CaseBrowserPanel::OpenCaseAction : " +
                                " no one picked up case"
                            );

                            // Fire case closed event.
                            fireCaseClosed
                            (
                                contextInfo.getCaseRoot(),
                                contextInfo.getCaseName()
                            );

                            // Nobody has picked up this case to manage
                            // so close it down.
                            caseServer.close();
                        }
                        else
                        {
                            // Mark the case as locked.
                            contextInfo.setCaseLocked(true);
                        }
                    }
                }
            }
            catch (FoamXIOError ioErr)
            {
                App.printMessage
                (
                    App.DEBUGLEVEL_VERBOSE,
                    "CaseBrowserPanel::OpenCaseAction : FoamXIOError"
                );

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

    private class OpenCasePostAction
        extends AbstractAction
    {
        OpenCasePostAction()
        {
            putValue
            (
                Action.SMALL_ICON, App.getResources().getIcon
                (
                    "CaseBrowser.OpenCaseImage"
                )
            );
            putValue(Action.NAME, "Postprocess Case");
            putValue(Action.SHORT_DESCRIPTION, "Postprocess Case");
            putValue(Action.LONG_DESCRIPTION, "Postprocess Case");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(CaseBrowserPanel.this);

            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode nodeItem =
                        (DefaultMutableTreeNode)tp.getLastPathComponent(); 

                    ContextInfo contextInfo = 
                        (ContextInfo)nodeItem.getUserObject();

                    // Check node type.
                    if (contextInfo.isCaseName())
                    {
                        // Get case browser reference from context info object.
                        ICaseBrowser caseBrowser = contextInfo.getCaseBrowser();
                        CaseDescriptor caseDescriptor =
                            contextInfo.getCaseDescriptor();

                        // Open the selected case.
                        ICasePostServerHolder caseHolder =
                            new ICasePostServerHolder();

                        int nProcs = caseDescriptor.nProcs;
                        if (nProcs != 1)
                        {
                            if (!openParallelDialog(nProcs))
                            {
                                nProcs = 1;
                            }
                        }

                        try
                        {
                            getCasePostServerReference
                            (
                                contextInfo.getCaseRoot(),
                                contextInfo.getCaseName(),
                                caseDescriptor,
                                nProcs,
                                caseBrowser,
                                caseHolder
                            );

                            //// Update status from C++ side.
                            //refreshCaseNode(contextInfo.getCaseRoot(), contextInfo.getCaseName());
                        }
                        catch (FoamXIOError ioErr)
                        {
                            // Mark the case as in error.
                            caseBrowser.caseIsInError(caseDescriptor);
                            refreshCaseNode(contextInfo.getCaseRoot(), contextInfo.getCaseName());

                            ioErr.errorMessage =
                                "The Open Case operation failed due to\n"
                              + ioErr.errorMessage;
                            throw ioErr;
                        }
                        catch (FoamXError fxErr)
                        {
                            // Mark the case as in error.
                            caseBrowser.caseIsInError(caseDescriptor);
                            refreshCaseNode(contextInfo.getCaseRoot(), contextInfo.getCaseName());

                            fxErr.errorMessage =
                                "The Open Case operation failed due to\n"
                              + fxErr.errorMessage;
                            throw fxErr;
                        }

                        ICasePostServer casePostServer = caseHolder.value;

                        // Fire case opened event.
                        fireCasePostOpened
                        (
                            contextInfo.getCaseRoot(),
                            contextInfo.getCaseName(),
                            casePostServer,
                            caseBrowser
                        );
                    }
                }
            }
            catch (FoamXIOError ioErr)
            {
                App.printMessage
                (
                    App.DEBUGLEVEL_VERBOSE,
                    "CaseBrowserPanel::OpenCasePostAction : FoamXIOError"
                );

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

    private class CreateCaseAction
    extends AbstractAction
    {
        CreateCaseAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("CaseBrowser.CreateCaseImage")
            );
            putValue(Action.NAME, "Create Case");
            putValue(Action.SHORT_DESCRIPTION, "Create Case");
            putValue(Action.LONG_DESCRIPTION, "Create Case");
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode nodeItem =
                        (DefaultMutableTreeNode)tp.getLastPathComponent();

                    ContextInfo contextInfo =
                        (ContextInfo)nodeItem.getUserObject();

                    // Check node type.
                    if (contextInfo.isHost() || contextInfo.isCaseRoot())
                    {
                        // Get case browser reference from context info object.
                        ICaseBrowser caseBrowser = contextInfo.getCaseBrowser();

                        // Get the available application classes
                        // for this browser.
                        ApplicationDescriptor[] appes =
                            caseBrowser.foamProperties().applicationes();

                        String[] caseRoots =
                            caseBrowser.foamProperties().rootDirectories();

                        // Create the create/import case dialog.
                        CreateCaseDialog dlg = new CreateCaseDialog
                        (
                            App.getRootFrame(),
                            appes,
                            caseRoots
                        );

                        dlg.setTitle("Create New Case");

                        if (contextInfo.isCaseRoot())
                        {
                            dlg.setCaseRoot(contextInfo.getCaseRoot());
                        }

                        // Show the dialog.
                        dlg.setVisible(true);

                        if (!dlg.wasCancelled())
                        {
                            // Show busy cursor.
                            BusyCursor busyCursor = new BusyCursor
                            (
                                CaseBrowserPanel.this
                            );

                            // Get the application class and case root and name.
                            String app = dlg.getCaseClass();
                            String caseRoot = dlg.getCaseRoot();
                            String caseName = dlg.getCaseName();

                            // Try and create a new case.
                            ICaseServerHolder caseHolder =
                                new ICaseServerHolder();

                            try
                            {
                                caseBrowser.newCase
                                (
                                    caseRoot,
                                    caseName,
                                    app
                                );

                                getCaseServerReference
                                (
                                    caseRoot,
                                    caseName,
                                    null,
                                    caseBrowser,
                                    caseHolder
                                );
                            }
                            catch (FoamXIOError ioErr)
                            {
                                //contextInfo.setCaseError(true);
                                ioErr.errorMessage =
                                    "The Create Case operation failed due to\n"
                                  + ioErr.errorMessage;
                                throw ioErr;
                            }
                            catch (FoamXError fxErr)
                            {
                                //contextInfo.setCaseError(true);
                                fxErr.errorMessage =
                                    "The Create Case operation failed due to\n"
                                  + fxErr.errorMessage;
                                throw fxErr;
                            }
                            ICaseServer caseServer = caseHolder.value;


                            // Fire case opened event.
                            fireCaseOpened
                            (
                                 caseRoot,
                                 caseName,
                                 caseServer,
                                 caseBrowser
                            );

                            if (!caseServer.managed())
                            {
                                // Fire case closed event.
                                fireCaseClosed(caseRoot, caseName);

                                // Nobody has picked up this case
                                // to manage so close it down.
                                caseServer.close();
                            }

                            // Refresh the case browser tree model.
                            refreshRoot(caseBrowser, caseRoot);
                        }
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

    private class DeleteCaseAction
    extends AbstractAction
    {
        DeleteCaseAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("CaseBrowser.DeleteCaseImage")
            );
            putValue(Action.NAME, "Delete Case");
            putValue(Action.SHORT_DESCRIPTION, "Delete Case");
            putValue(Action.LONG_DESCRIPTION, "Delete Case");
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode nodeItem = 
                        (DefaultMutableTreeNode)tp.getLastPathComponent();
                    ContextInfo contextInfo =
                        (ContextInfo)nodeItem.getUserObject();

                    // Check node type.
                    if (contextInfo.isCaseName())
                    {
                        // Get case browser reference from context info object.
                        ICaseBrowser caseBrowser = contextInfo.getCaseBrowser();
                        CaseDescriptor caseDescriptor =
                            contextInfo.getCaseDescriptor();

                        // Final check to see if the case is locked.
                        if (caseBrowser.caseLocked(caseDescriptor))
                        {
                            // Issue warning.
                            JOptionPane.showMessageDialog
                            (
                                App.getRootFrame(),
                                "This case is locked."
                              + "\n(it is probably already being preprocessed)",
                                "Case Locked Error",
                                JOptionPane.ERROR_MESSAGE
                            );

                            // Mark as locked.
                            contextInfo.setCaseLocked(true);
                            caseBrowser = null;
                            return;
                        }

                        String caseName = contextInfo.getCaseName();
                        String caseRoot = contextInfo.getCaseRoot();

                        // Prompt the user to see if the case files are
                        // to be deleted.
                        int ret = JOptionPane.showConfirmDialog
                        (
                            App.getRootFrame(),
                            "Delete case " + caseName + "?"
                          + "\n(deletes case directory)",
                            "Confirm Delete Case Operation",
                            JOptionPane.YES_NO_OPTION,
                            JOptionPane.QUESTION_MESSAGE
                        );

                        // See if the user has changed his mind. Bloody users.
                        if (ret == JOptionPane.OK_OPTION)
                        {
                            // Show busy cursor.
                            BusyCursor busyCursor = new BusyCursor
                            (
                                CaseBrowserPanel.this
                            );

                            // Try and delete the case.
                            try
                            {
                                caseBrowser.deleteCase
                                (
                                    contextInfo.getCaseDescriptor()
                                );
                            }
                            catch (FoamXIOError ioErr)
                            {
//                                caseBrowser.caseIsInError(caseDescriptor);
//                                contextInfo.setCaseError(true);

                                ioErr.errorMessage =
                                    "The Delete Case operation failed due to\n"
                                  + ioErr.errorMessage;
                                throw ioErr;
                            }
                            catch (FoamXError fxErr)
                            {
//                                caseBrowser.caseIsInError(caseDescriptor);
//                                contextInfo.setCaseError(true);

                                fxErr.errorMessage =
                                    "The Delete Case operation failed due to\n"
                                  + fxErr.errorMessage;
                                throw fxErr;
                            }

                            // Fire case deleted event.
                            fireCaseDeleted
                            (
                                caseRoot,
                                caseName
                            );


//                            // Delete case operation has been successful
//                            //so delete the case node.
//                            model_.removeNodeFromParent(nodeItem);
//                            // Refresh the case browser tree model.
//                            tree_.updateUI();

                            // Refresh the case browser tree model for this
                            // caseRoot.
                            refreshRoot(caseBrowser, caseRoot);
                        }
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

    private class CloneCaseAction
    extends AbstractAction
    {
        CloneCaseAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("CaseBrowser.CloneCaseImage")
            );
            putValue(Action.NAME, "Clone Case");
            putValue(Action.SHORT_DESCRIPTION, "Clone Case");
            putValue(Action.LONG_DESCRIPTION, "Clone Case");
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode nodeItem = 
                        (DefaultMutableTreeNode)tp.getLastPathComponent();
                    ContextInfo contextInfo =
                        (ContextInfo)nodeItem.getUserObject();

                    // Check node type.
                    if (contextInfo.getType() == ContextInfo.CASENAME)
                    {
                        // Get case browser reference from context info object.
                        ICaseBrowser caseBrowser = contextInfo.getCaseBrowser();
                        CaseDescriptor caseDescriptor =
                            contextInfo.getCaseDescriptor();

                        // Final check to see if the case is locked.
                        if (caseBrowser.caseLocked(caseDescriptor))
                        {
                            // Issue warning.
                            JOptionPane.showMessageDialog
                            (
                                App.getRootFrame(),
                                "This case is locked."
                              + "\n(it is probably already being preprocessed)",
                                "Case Locked Error",
                                JOptionPane.ERROR_MESSAGE
                            );

                            // Mark as locked.
                            contextInfo.setCaseLocked(true);
                            caseBrowser = null;
                            return;
                        }

                        // Get the available application classes
                        // for this browser.
                        ApplicationDescriptor[] appes =
                            caseBrowser.foamProperties().applicationes();

                        IDictionaryEntryHolder dictHolder =
                            new IDictionaryEntryHolder();
                        caseBrowser.foamProperties().getUtilityControlDict
                        (
                            "foamCloneCase",
                            contextInfo.getCaseRoot(),
                            contextInfo.getCaseName(),
                            dictHolder
                        );

                        if (dictHolder.value == null)
                        {
                            throw new FoamXException
                            (
                                "Could not find dictionary for foamCloneCase"
                            );
                        }

                        // Get dictionary
                        DictionaryEntryCache controlDict =
                            DictionaryEntryCache.New
                            (
                                dictHolder.value
                            );

                        // Set default root
                        DictionaryEntryCache rootEntry =
                            controlDict.getSubEntry
                            (
                                "root"
                            );
                        rootEntry.updateValue(contextInfo.getCaseRoot());

                        // Convert available appclasses into stringlist
                        String[] appNames = new String[appes.length];
                        for (int i = 0; i < appes.length; i++)
                        {
                            appNames[i] = appes[i].name;
                        }
                        DictionaryEntryCache appEntry =
                            controlDict.getSubEntry
                            (
                                "application"
                            );

                        // Update 'valueList' of appEntry and write
                        // to case browser.
                        appEntry.getTypeDescriptor().setValueList
                        (
                            appNames
                        );
                        appEntry.getTypeDescriptor().updateTypeDescriptor
                        (
                            true
                        );
                        appEntry.updateValue(caseDescriptor.app);

                        CompoundEditor compoundEditor =
                            new CompoundEditor
                            (
                                App.getRootFrame(),
                                controlDict.getDictEntry()
                            );
                        compoundEditor.setTitle("Clone Case");
                        compoundEditor.setVisible(true);

                        // Show busy cursor.
                        BusyCursor busyCursor = new BusyCursor
                        (
                            CaseBrowserPanel.this
                        );

                        DictionaryEntryCache rootDirEntry =
                        controlDict.getSubEntry
                        (
                            "root"
                        );
                        String caseRoot = rootDirEntry.toString();

                        DictionaryEntryCache caseNameEntry =
                        controlDict.getSubEntry
                        (
                            "case"
                        );
                        String caseName = caseNameEntry.toString();

                        DictionaryEntryCache timesEntry =
                        controlDict.getSubEntry
                        (
                            "times"
                        );
                        String times = timesEntry.toString();

                        DictionaryEntryCache appNameEntry =
                        controlDict.getSubEntry
                        (
                            "application"
                        );
                        String app = appNameEntry.toString();


                        // Pop up confirmation box
                        if
                        (
                            JOptionPane.showConfirmDialog
                            (
                                App.getRootFrame(),
                                "Clone from "
                              + caseDescriptor.caseName
                              + " to case "
                              + caseName + " with time directories "
                              + times + "?",
                                "Confirm Clone Case Operation",
                                JOptionPane.YES_NO_OPTION,
                                JOptionPane.QUESTION_MESSAGE
                            )
                            ==
                            JOptionPane.OK_OPTION
                        )
                        {
                            // Try and clone the case.
                            try
                            {
                                caseBrowser.cloneCase
                                (
                                    caseDescriptor,
                                    caseRoot,
                                    caseName,
                                    app,
                                    times
                                );
                            }
                            catch (FoamXIOError ioErr)
                            {
                                ioErr.errorMessage =
                                    "The Clone Case operation failed due to\n"
                                  + ioErr.errorMessage;
                                throw ioErr;
                            }
                            catch (FoamXError fxErr)
                            {
                                fxErr.errorMessage =
                                    "The Clone Case operation failed due to\n"
                                  + fxErr.errorMessage;
                                throw fxErr;
                            }

                            // Refresh the case browser for this caseRoot.
                            refreshRoot
                            (
                                caseBrowser,
                                contextInfo.getCaseRoot()
                            );
                        }
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

    private class UnlockCaseAction
    extends AbstractAction
    {
        UnlockCaseAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("CaseBrowser.UnlockCaseImage")
            );
            putValue(Action.NAME, "Unlock Case");
            putValue(Action.SHORT_DESCRIPTION, "Unlock Case");
            putValue(Action.LONG_DESCRIPTION, "Unlock Case");
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode nodeItem =
                        (DefaultMutableTreeNode)tp.getLastPathComponent();
                    ContextInfo contextInfo =
                        (ContextInfo)nodeItem.getUserObject();

                    // Check node type.
                    if (contextInfo.getType() == ContextInfo.CASENAME)
                    {
                        // Get case browser reference from context info object.
                        ICaseBrowser caseBrowser = contextInfo.getCaseBrowser();

                        // Final check to see if the case is really locked.
                        if
                        (
                            caseBrowser.caseLocked
                            (
                                contextInfo.getCaseDescriptor()
                            )
                        )
                        {
                            // Issue warning.
                            if
                            (
                                JOptionPane.showConfirmDialog
                                (
                                    App.getRootFrame(),
                                    "Case is currently locked.\n"
                                    + "Do you wish to unlock it?",
                                    "Confirm Unlock Case Operation",
                                    JOptionPane.YES_NO_OPTION,
                                    JOptionPane.QUESTION_MESSAGE
                                )
                                ==
                                JOptionPane.OK_OPTION
                            )
                            {
                                // Unlock the case.
                                caseBrowser.unlockCaseDescriptor
                                (
                                    contextInfo.getCaseDescriptor()
                                );
//                                // Mark as unlocked.
//                                contextInfo.setCaseLocked(false);
//                                contextInfo.setCaseError(false);

                                // Update case browser.
                                refreshCaseNode
                                (
                                    contextInfo.getCaseRoot(),
                                    contextInfo.getCaseName()
                                );
                            }
                        }
                        else
                        {
                            // Case not locked on disk but might be in
                            // caseBrowser.
                            caseBrowser.unlockCaseDescriptor
                            (
                                contextInfo.getCaseDescriptor()
                            );

                            // Update case browser.
                            refreshCaseNode
                            (
                                contextInfo.getCaseRoot(),
                                contextInfo.getCaseName()
                            );
//                            // Mark as unlocked.
//                            contextInfo.setCaseLocked(false);
//                            contextInfo.setCaseError(false);
                        }
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

    private class RefreshAction
    extends AbstractAction
    {
        RefreshAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("CaseBrowser.RefreshBrowserImage")
            );
            putValue(Action.NAME, "Refresh Case Browser");
            putValue(Action.SHORT_DESCRIPTION, "Refresh Case Browser");
            putValue(Action.LONG_DESCRIPTION, "Refresh Case Browser");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(App.getRootFrame());

            try
            {
                // Refresh the case browser tree model.
                if
                (
                    JOptionPane.showConfirmDialog
                    (
                        App.getRootFrame(),
                        "Ok to refresh hosts and caselists?"
                        + "\n(This is done in the background)",
                        "Refresh Hosts and Cases",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.QUESTION_MESSAGE
                    ) ==
                    JOptionPane.OK_OPTION
                )
                {
                    // Do the checking in separate thread.
                    Thread refreshThread = new Thread()
                    {
                        public void run()
                        {
                            refreshAll();
                        }
                    };
                    refreshThread.start();
                }
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }


    //--------------------------------------------------------------------------
    //---- Action Classes for Selection Mode
    //--------------------------------------------------------------------------

    private class SelectCaseAction
        extends AbstractAction
    {
        SelectCaseAction()
        {
            putValue
            (
                Action.SMALL_ICON, App.getResources().getIcon
                (
                    "CaseBrowser.OpenCaseImage"
                )
            );
            putValue(Action.NAME, "Select Case");
            putValue(Action.SHORT_DESCRIPTION, "Select Case");
            putValue(Action.LONG_DESCRIPTION, "Select Case");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(CaseBrowserPanel.this);

            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    DefaultMutableTreeNode aaa =
                        (DefaultMutableTreeNode)tp.getPathComponent(0);

                    ContextInfo bbb = 
                        (ContextInfo)aaa.getUserObject();

                    // Get node info.
                    DefaultMutableTreeNode nodeItem =
                        (DefaultMutableTreeNode)tp.getLastPathComponent(); 

                    ContextInfo contextInfo = 
                        (ContextInfo)nodeItem.getUserObject();

                    // Check node type.
                    if (contextInfo.isCaseName())
                    {
                        // Get case browser reference from context info object.
                        ICaseBrowser caseBrowser =
                            contextInfo.getCaseBrowser();

                        StringHolder holder = new StringHolder();
                        caseBrowser.getHostName(holder);

                        // Fire case selected event.
                        fireCaseSelected
                        (
                            holder.value,
                            contextInfo.getCaseRoot(),
                            contextInfo.getCaseName(),
                            caseBrowser
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

    private class SelectRootAction
        extends AbstractAction
    {
        SelectRootAction()
        {
            putValue
            (
                Action.SMALL_ICON, App.getResources().getIcon
                (
                    "CaseBrowser.OpenCaseImage"
                )
            );
            putValue(Action.NAME, "Select Root");
            putValue(Action.SHORT_DESCRIPTION, "Select Root");
            putValue(Action.LONG_DESCRIPTION, "Select Root");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(CaseBrowserPanel.this);

            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode nodeItem =
                        (DefaultMutableTreeNode)tp.getLastPathComponent(); 

                    ContextInfo contextInfo = 
                        (ContextInfo)nodeItem.getUserObject();

                    // Check node type.
                    if (contextInfo.isCaseRoot())
                    {
                        // Get case browser reference from context info object.
                        ICaseBrowser caseBrowser =
                            contextInfo.getCaseBrowser();

                        StringHolder holder = new StringHolder();
                        caseBrowser.getHostName(holder);

                        // Fire root selected event.
                        fireRootSelected
                        (
                            holder.value,
                            contextInfo.getCaseRoot(),
                            caseBrowser
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

    private class SelectHostAction
        extends AbstractAction
    {
        SelectHostAction()
        {
            putValue
            (
                Action.SMALL_ICON, App.getResources().getIcon
                (
                    "CaseBrowser.OpenCaseImage"
                )
            );
            putValue(Action.NAME, "Select Host");
            putValue(Action.SHORT_DESCRIPTION, "Select Host");
            putValue(Action.LONG_DESCRIPTION, "Select Host");
        }

        public void actionPerformed(ActionEvent e)
        {
            // Show busy cursor.
            BusyCursor busyCursor = new BusyCursor(CaseBrowserPanel.this);

            try
            {
                // Get currently selected node.
                TreePath tp = tree_.getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode nodeItem =
                        (DefaultMutableTreeNode)tp.getLastPathComponent(); 

                    ContextInfo contextInfo = 
                        (ContextInfo)nodeItem.getUserObject();

                    // Check node type.
                    if (contextInfo.isHost())
                    {
                        // Fire case opened event.
                        fireHostSelected(contextInfo.getText());
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
}
