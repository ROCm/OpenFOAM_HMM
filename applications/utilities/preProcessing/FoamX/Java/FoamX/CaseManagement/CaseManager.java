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

import java.util.Hashtable;
import java.util.Enumeration;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

import FoamX.App;
import FoamX.TaskManagement.Task;
import FoamX.ToolbarManagement.Toolbar;
import FoamX.Util.BusyCursor;
import FoamX.Exceptions.FoamXException;
import FoamX.Exceptions.CaseBrowserIOException;

import FoamXServer.*;
import FoamXServer.CaseServer.ICaseServer;
import FoamXServer.CasePostServer.ICasePostServer;

public class CaseManager
    extends javax.swing.JPanel
    implements CaseStatusListener, ChangeListener
{
    //--------------------------------------------------------------------------

    // Case browser panel object.
    protected CaseBrowserPanel caseBrowserPanel_; 
    // Map of CasePanel objects.
    protected Hashtable caseMap_;
    // Map of CasePostPanel objects.
    protected Hashtable casePostMap_;
    // Objects interested in case manager events.
    protected EventListenerList listenerList_;
    protected Toolbar toolBar_;
    //protected SaveCurrentCaseAction saveCurrentCaseAction_;
    //protected CloseCurrentCaseAction closeCurrentCaseAction_;
    //protected SaveAllCasesAction saveAllCasesAction_;
    //protected CloseAllCasesAction closeAllCasesAction_;

    //--------------------------------------------------------------------------
    /** CaseManager constructor. */
    public CaseManager()
    {
        try
        {
            // Initialise the GUI components.
            initComponents();

            // Create case map.
            caseMap_ = new Hashtable();
            casePostMap_ = new Hashtable();

            // Create listener list.
            listenerList_ = new EventListenerList();

            // Initialise actions and toolbar.
            //saveCurrentCaseAction_  = new SaveCurrentCaseAction();
            //closeCurrentCaseAction_ = new CloseCurrentCaseAction();
            //saveAllCasesAction_     = new SaveAllCasesAction();
            //closeAllCasesAction_    = new CloseAllCasesAction();

            //saveCurrentCaseAction_.setEnabled(false);
            //closeCurrentCaseAction_.setEnabled(false);
            //saveAllCasesAction_.setEnabled(false);
            //closeAllCasesAction_.setEnabled(false);


            // Add case browser tab. Have it start a new model.
            caseBrowserPanel_ =
                new CaseBrowserPanel(null, CaseBrowserPanel.OPEN_CASE_MODE);
            tabbedPanel_.addTab("Case Browser", caseBrowserPanel_);

            // Initialise toolbar with my actions
            //toolBar_ =
            //    App.getToolbarManager().addToolbar
            //    (
            //        "CaseManager",
            //        "Case Manager"
            //    );
            //addToolbarButton(saveCurrentCaseAction_);
            //addToolbarButton(saveAllCasesAction_);
            //toolBar_.addSeparator();
            //addToolbarButton(closeCurrentCaseAction_);
            //addToolbarButton(closeAllCasesAction_);

            // Register this object as a tab selection changed listener.
            tabbedPanel_.getModel().addChangeListener(this);

            // Register this object as a CaseStatusListener.
            caseBrowserPanel_.addCaseStatusListener(this);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Close all open cases. */
    public void closeAllCases()
    {
        try
        {
            // Create a task object to keep track of progress.
            Task task = App.getTaskManager().createTask
                (
                    "Closing All Cases",
                    0,
                    caseMap_.size()
                );
            int caseNo = 0;

            // Close all cases.
            Enumeration iter = caseMap_.elements();
            while (iter.hasMoreElements())
            {
                // Get next case.
                CasePanel casePanel = (CasePanel)iter.nextElement();

                // Give option of saving if modified
                casePanel.checkAndSave
                (
                    false,
                    "Case has been changed. Save before closing?",
                    "Save Case"
                );

                // Do exactly if 'close case' on casePanel was pressed
                caseClosed
                (
                    new CaseStatusEvent
                    (
                        this,
                        casePanel.getCaseRoot(),
                        casePanel.getCaseName()
                    )
                );

                // Update the progress info.
                task.setProgress(caseNo++);
            }


            //// Close all post processing in a similar way
            //task = App.getTaskManager().createTask
            //    (
            //        "Closing All Postprocessing",
            //        0,
            //        caseMap_.size()
            //    );
            //caseNo = 0;
            //// Close all cases.
            //iter = casePostMap_.elements();
            //while (iter.hasMoreElements())
            //{
            //    // Get next case.
            //    CasePostPanel casePostPanel = (CasePostPanel)iter.nextElement();
            //
            //    // Do exactly if 'close case' on casePostPanel was pressed
            //    casePostClosed
            //    (
            //        new CaseStatusEvent
            //        (
            //            this,
            //            casePanel.getCaseRoot(),
            //            casePanel.getCaseName()
            //        )
            //    );
            //
            //    // Update the progress info.
            //    task.setProgress(caseNo++);
            //}


            // Refresh the case browser panel
            caseBrowserPanel_.refreshCaseBrowser(false);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
        finally
        {
            App.getTaskManager().endTask();
        }
    }

    //--------------------------------------------------------------------------
    /** Save open case. Ask first. */
    public void saveCaseNice(String caseRoot, String caseName)
    {
        CasePanel casePanel = getCasePanel(caseRoot, caseName);

        if (casePanel != null)
        {
            // Give option of saving if modified
            casePanel.checkAndSave
            (
                false,
                "Case " + caseName
              + " has been changed. Save before invoking utility?",
                "Save Case"
            );

        }
    }


    //--------------------------------------------------------------------------
    /** Shut down the CaseManager object and close all open cases. */
    public void shutdown()
    {
        try
        {
            // Close all cases.
            closeAllCases();

            // Unregister as a CaseStatusListener.
            caseBrowserPanel_.removeCaseStatusListener(this);

            // Shutdown the case browser.
            caseBrowserPanel_.shutdown();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public CaseBrowserPanel getCaseBrowser()
    {
        return caseBrowserPanel_;
    }

    //--------------------------------------------------------------------------

    public CasePanel getCasePanel(String caseRoot, String caseName)
    {
        CasePanel casePanel = null;

        String rootAndCase = caseRoot + "/" + caseName;

        if (caseMap_.containsKey(rootAndCase))
        {
            // Create a new case panel object to manage the case.
            casePanel = (CasePanel)caseMap_.get(rootAndCase);
        }

        return casePanel;
    }

    //--------------------------------------------------------------------------

    protected void addToolbarButton(AbstractAction action)
    {
        // Add new button to the toolbar using the action object.
        javax.swing.JButton button = toolBar_.add(action);

        // Set the tooltip text.
        button.setToolTipText((String)action.getValue
            (
                Action.SHORT_DESCRIPTION)
            );
        button.setText("");
        button.setFont(toolBar_.getFont());      // Use same font as toolbar.
    }

    //--------------------------------------------------------------------------
    //---- CaseManagerListener Interface
    //--------------------------------------------------------------------------

    public void addCaseManagerListener(CaseManagerListener l)
    {
        listenerList_.add(CaseManagerListener.class, l);
    }

    //--------------------------------------------------------------------------

    public void removeCaseManagerListener(CaseManagerListener l)
    {
        listenerList_.remove(CaseManagerListener.class, l);
    }

    //--------------------------------------------------------------------------

    protected void fireOpenCase
    (
        String caseRoot,
        String caseName,
        ICaseServer caseServer
    )
    {
        // Create event object.
        CaseManagerEvent evt = new CaseManagerEvent
            (
                this,
                caseRoot,
                caseName,
                caseServer
            );

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseManagerListener.class)
            {
                ((CaseManagerListener)listeners[i+1]).caseOpened(evt);
            }
        }
    }


    //--------------------------------------------------------------------------

    protected void fireOpenPostCase
    (
        String caseRoot,
        String caseName,
        ICasePostServer casePostServer
    )
    {
        // Create event object.
        CaseManagerEvent evt = new CaseManagerEvent
            (
                this,
                caseRoot,
                caseName,
                casePostServer
            );

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseManagerListener.class)
            {
                ((CaseManagerListener)listeners[i+1]).casePostOpened(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected void fireCloseCase(String caseRoot, String caseName)
    {
        // Create event object.
        CaseManagerEvent evt = new CaseManagerEvent(this, caseRoot, caseName);

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseManagerListener.class)
            {
                ((CaseManagerListener)listeners[i+1]).caseClosed(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected void fireClosePostCase(String caseRoot, String caseName)
    {
        // Create event object.
        CaseManagerEvent evt = new CaseManagerEvent(this, caseRoot, caseName);

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseManagerListener.class)
            {
                ((CaseManagerListener)listeners[i+1]).casePostClosed(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected void fireDeleteCase(String caseRoot, String caseName)
    {
        // Create event object.
        CaseManagerEvent evt = new CaseManagerEvent(this, caseRoot, caseName);

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseManagerListener.class)
            {
                ((CaseManagerListener)listeners[i+1]).caseDeleted(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected void fireCasePanelSelected(String caseRoot, String caseName)
    {
        // Create event object.
        CaseManagerEvent evt = new CaseManagerEvent(this, caseRoot, caseName);

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseManagerListener.class)
            {
                ((CaseManagerListener)listeners[i+1]).casePanelSelected(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected void fireCaseBrowserSelected()
    {
        // Create event object.
        CaseManagerEvent evt = new CaseManagerEvent(this);

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == CaseManagerListener.class)
            {
                ((CaseManagerListener)listeners[i+1]).caseBrowserSelected(evt);
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
        tabbedPanel_ = new javax.swing.JTabbedPane();
        
        setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        tabbedPanel_.setTabPlacement(javax.swing.JTabbedPane.BOTTOM);
        tabbedPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        add(tabbedPanel_, gridBagConstraints1);
        
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JTabbedPane tabbedPanel_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //---- CaseStatusListener Interface
    //--------------------------------------------------------------------------

    public void caseOpened(CaseStatusEvent evt)
    {
        try
        {
            // Make sure the event object has a case server reference.
            if (evt.caseServer() != null)
            {
                // See if we already have this case.
                String key = evt.caseRoot() + "/" + evt.caseName();

                if (!caseMap_.containsKey(key))
                {
                    // Create a new case panel object to manage the case.
                    // Note: will throw exception if case in error.
                    CasePanel casePanel = new CasePanel
                    (
                        evt.caseBrowser(),
                        evt.caseServer()
                    );

                    // Register this object as a CaseStatusListener.
                    casePanel.addCaseStatusListener(this);


                    // Add new tab.
                    tabbedPanel_.addTab(evt.caseName(), casePanel);
                    caseMap_.put(key, casePanel);

                    // Pop the new tab to the front.
                    tabbedPanel_.setSelectedComponent(casePanel);

                    // Inform case manager listeners.
                    fireOpenCase
                    (
                        evt.caseRoot(),
                        evt.caseName(),
                        evt.caseServer()
                    );
                }
            }
        }
        catch (FoamXIOError ioErr)
        {
            // io error. Add caseBrowser information
            if (evt.caseBrowser() != null)
            {
                // Add caseBrowser info
                App.handleAllExceptions
                (
                    new CaseBrowserIOException
                    (
                        ioErr,
                        evt.caseBrowser()
                    )
                );
            }
            else
            {
                App.handleAllExceptions(ioErr);
            }
        }
        catch (Exception ex)
        {
//            caseBrowserPanel_.refreshCaseNode
//            (
//                evt.caseRoot(),
//                evt.caseName()
//            );
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void casePostOpened(CaseStatusEvent evt)
    {
        try
        {
            // Make sure the event object has a case server reference.
            if (evt.casePostServer() == null)
            {
                return;
            }

            // See if we already have this case.
            String key = evt.caseRoot() + "/" + evt.caseName();

            if (casePostMap_.containsKey(key))
            {
                System.out.println("Already postprocessing case " + key);
            }
            else
            {
                System.out.println("Starting postprocessing case " + key);
            }

            //PostApp pp =
            //    new PostApp
            //    (
            //        evt.casePostServer(),
            //        evt.caseRoot(),
            //        evt.caseName()
            //    );
            ////MainFrame mf = new MainFrame(pp, 300, 600);
            //
            //// Register this object as a CaseStatusListener.
            //pp.getPostWindow().addCaseStatusListener(this);
            //
            //casePostMap_.put(key, pp);

            // Inform case manager listeners.
            fireOpenPostCase
            (
                evt.caseRoot(),
                evt.caseName(),
                evt.casePostServer()
            );
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void caseClosed(CaseStatusEvent evt)
    {
        try
        {
            // See if we have this case.
            String key = evt.caseRoot() + "/" + evt.caseName();

            if (caseMap_.containsKey(key))
            {
                // Get CasePanel object for this case.
                CasePanel casePanel = (CasePanel)caseMap_.get(key);

                // Unregister as a CaseStatusListener.
                casePanel.removeCaseStatusListener(this);

                // Close the case.
                casePanel.closeCase(CasePanel.KILL_SERVER);

                // Remove tab.
                tabbedPanel_.remove(casePanel);
                caseMap_.remove(key);

                // Pop the case browser tab to the front.
                tabbedPanel_.setSelectedComponent(caseBrowserPanel_);

                // Refresh the case browser panel
                caseBrowserPanel_.refreshCaseNode
                (
                    evt.caseRoot(),
                    evt.caseName()
                );

                // Inform case manager listeners.
                fireCloseCase(evt.caseRoot(), evt.caseName());
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void casePostClosed(CaseStatusEvent evt)
    {
        try
        {
            //// See if we have this case.
            //String key = evt.caseRoot() + "/" + evt.caseName();
            //
            //if (casePostMap_.containsKey(key))
            //{
            //    System.out.println
            //    (
            //        "CaseManager:casePostClosed for " + evt.caseName()
            //    );
            //
            //    PostApp pp = (PostApp)casePostMap_.get(key);
            //
            //    PostWindow ps = pp.getPostWindow();
            //
            //    // Unregister as a CaseStatusListener.
            //    ps.removeCaseStatusListener(this);
            //
            //    pp.destroy();
            //
            //    casePostMap_.remove(key);
            //
            //    // Inform case manager listeners.
            //    fireClosePostCase(evt.caseRoot(), evt.caseName());
            //}
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void caseDeleted(CaseStatusEvent evt)
    {
        try
        {
            // See if we have this case.
            String key = evt.caseRoot() + "/" + evt.caseName();
            if (caseMap_.containsKey(key))
            {
                // Get CasePanel object for this case.
                CasePanel casePanel = (CasePanel)caseMap_.get(key);

                // Unregister as a CaseStatusListener.
                casePanel.removeCaseStatusListener(this);

                // Close the case.
                casePanel.closeCase(CasePanel.KILL_SERVER);

                // Remove tab.
                tabbedPanel_.remove(casePanel);
                caseMap_.remove(key);

                // Pop the case browser tab to the front.
                tabbedPanel_.setSelectedComponent(caseBrowserPanel_);

                // Refresh the case browser.
                caseBrowserPanel_.deleteCaseNode
                (
                    evt.caseRoot(),
                    evt.caseName()
                );

                // Inform case manager listeners.
                fireDeleteCase(evt.caseRoot(), evt.caseName());
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //---- ChangeEventListener Interface
    //--------------------------------------------------------------------------

    public void stateChanged(javax.swing.event.ChangeEvent evt)
    {
        try
        {
            // See which tab has been selected;
            Component comp = tabbedPanel_.getSelectedComponent();
            if (comp instanceof CaseBrowserPanel)
            {
                // Tell the case browser tab that it's selected.
                caseBrowserPanel_.tabSelected();

                // Deselect all case tabs.
                Component[] tabs = tabbedPanel_.getComponents();
                for (int i = 0; i <tabs.length; i++)
                {
                    if (tabs[i] instanceof CasePanel)
                    {
                        CasePanel cp = (CasePanel)tabs[i];
                        cp.tabDeselected();
                    }
                }

                // Fire a tab selection event to interested listeners.
                fireCaseBrowserSelected();
            }
            else if (comp instanceof CasePanel)
            {
                // Deselect all case tabs and the case browser tab.
                caseBrowserPanel_.tabDeselected();
                Component[] tabs = tabbedPanel_.getComponents();
                for (int i = 0; i <tabs.length; i++)
                {
                    if (tabs[i] instanceof CasePanel)
                    {
                        CasePanel cp = (CasePanel)tabs[i];
                        cp.tabDeselected();
                    }
                }

                // Tell the case panel that it's selected.
                CasePanel cp = (CasePanel)comp;
                cp.tabSelected();

                // Fire a tab selection event to interested listeners.
                fireCasePanelSelected(cp.getCaseRoot(), cp.getCaseName());
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //---- Action Classes
    //--------------------------------------------------------------------------

//    private class SaveCurrentCaseAction
//    extends AbstractAction
//    {
//        SaveCurrentCaseAction()
//        {
//            putValue
//            (
//                Action.SMALL_ICON,
//                App.getResources().getIcon
//                    (
//                        "CaseManager.SaveCurrentCaseImage"
//                    )
//            );
//            putValue(Action.NAME, "SaveCurrentCase");
//            putValue(Action.SHORT_DESCRIPTION, "Save Current Case");
//            putValue(Action.LONG_DESCRIPTION, "Save Current Case");
//        }
//
//        public void actionPerformed(ActionEvent evt)
//        {
//            // Show busy cursor.
//            BusyCursor cursor = new BusyCursor(App.getRootFrame());
//
//            try
//            {
//                Component comp = tabbedPanel_.getSelectedComponent();
//                if (comp instanceof CasePanel)
//                {
//                    CasePanel cp = (CasePanel)comp;
//                    cp.saveCase();
//                }
//            }
//            catch (Exception ex)
//            {
//                App.handleAllExceptions(ex);
//            }
//        }
//    }

    //--------------------------------------------------------------------------

//    private class SaveAllCasesAction
//    extends AbstractAction
//    {
//        SaveAllCasesAction()
//        {
//            putValue
//            (
//                Action.SMALL_ICON,
//                App.getResources().getIcon("CaseManager.SaveAllCasesImage")
//            );
//            putValue(Action.NAME, "SaveAllCases");
//            putValue(Action.SHORT_DESCRIPTION, "Save All Cases");
//            putValue(Action.LONG_DESCRIPTION, "Save All Cases");
//        }
//
//        public void actionPerformed(ActionEvent evt)
//        {
//            // Show busy cursor.
//            BusyCursor cursor = new BusyCursor(App.getRootFrame());
//
//            try
//            {
//                Component[] tabs = tabbedPanel_.getComponents();
//                for (int i = 0; i <tabs.length; i++)
//                {
//                    if (tabs[i] instanceof CasePanel)
//                    {
//                        CasePanel cp = (CasePanel)tabs[i];
//                        cp.saveCase();
//                    }
//                }
//            }
//            catch (Exception ex)
//            {
//                App.handleAllExceptions(ex);
//            }
//        }
//    }

    //--------------------------------------------------------------------------

//    private class CloseCurrentCaseAction
//    extends AbstractAction
//    {
//        CloseCurrentCaseAction()
//        {
//            putValue
//            (
//                Action.SMALL_ICON,
//                App.getResources().getIcon("CaseManager.CloseCurrentCaseImage")
//            );
//            putValue(Action.NAME, "CloseCurrentCase");
//            putValue(Action.SHORT_DESCRIPTION, "Close Current Case");
//            putValue(Action.LONG_DESCRIPTION, "Close Current Case");
//        }
//
//        public void actionPerformed(ActionEvent evt)
//        {
//            // Show busy cursor.
//            BusyCursor cursor = new BusyCursor(App.getRootFrame());
//
//            try
//            {
//                Component comp = tabbedPanel_.getSelectedComponent();
//                if (comp instanceof CasePanel)
//                {
//                    CasePanel cp = (CasePanel)comp;
//                    boolean doClose = cp.checkAndSave
//                    (
//                        true,
//                        "Case has been changed. Save before closing?",
//                        "Save Case"
//                    );
//                    if (doClose)
//                    {
//                        // Do exactly if 'close case' on casePanel was pressed
//                        caseClosed
//                        (
//                            new CaseStatusEvent
//                            (
//                                this,
//                                cp.getCaseRoot(),
//                                cp.getCaseName()
//                            )
//                        );
//                        
//                    }
//                }
//            }
//            catch (Exception ex)
//            {
//                App.handleAllExceptions(ex);
//            }
//        }
//    }

    //--------------------------------------------------------------------------

//    private class CloseAllCasesAction
//    extends AbstractAction
//    {
//        CloseAllCasesAction()
//        {
//            putValue
//            (
//                Action.SMALL_ICON,  
//                App.getResources().getIcon("CaseManager.CloseAllCasesImage")
//            );
//            putValue(Action.NAME, "CloseAllCases");
//            putValue(Action.SHORT_DESCRIPTION, "Close All Cases");
//            putValue(Action.LONG_DESCRIPTION, "Close All Cases");
//        }
//
//        public void actionPerformed(ActionEvent evt)
//        {
//            // Show busy cursor.
//            BusyCursor cursor = new BusyCursor(App.getRootFrame());
//
//            try
//            {
//                closeAllCases();
//            }
//            catch (Exception ex)
//            {
//                App.handleAllExceptions(ex);
//            }
//        }
//    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}




