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
import java.util.Properties;
import java.util.Enumeration;
import javax.swing.*;
import javax.swing.tree.*;

import org.omg.CosNaming.*;

import FoamX.App;
import FoamX.Exceptions.FoamXException;

import FoamXServer.*;
import FoamXServer.CaseServer.*;

import FoamXServer.CaseBrowser.ICaseBrowser;
import FoamXServer.CaseBrowser.ICaseBrowserHolder;
import FoamXServer.CaseDescriptor;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;
import FoamXServer.HostBrowser.IHostBrowser;

public class CaseBrowserModel
    extends DefaultTreeModel
{
    //--------------------------------------------------------------------------

    protected IHostBrowser hostBrowser_;
    protected DefaultMutableTreeNode root_;
    protected Hashtable hostMap_;       // from host to treeNode
    protected Hashtable rootDirMap_;    // from root to treeNode
    protected Hashtable keyMap_;        // from key (root/case) to treeNode
    protected ICaseBrowser caseBrowserCopy_;

    //--------------------------------------------------------------------------
    /** CaseManagerModel constructor. */
    public CaseBrowserModel()
    {
        super(new DefaultMutableTreeNode());

        try
        {
            // Get reference to root node.
            root_ = (DefaultMutableTreeNode)getRoot();

            // Create host map.
            hostMap_ = new Hashtable(10);
            rootDirMap_ = new Hashtable(10);
            keyMap_ = new Hashtable(10);

            // Initialise the model.
            // Show only the Hosts node unless the model is refreshed.
            hostBrowser_ = null;
            root_.setUserObject
            (
                new ContextInfo("Hosts", ContextInfo.INVALID_REF)
            );
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Initialise the (licenseHosts) tree model. */
    protected void initialiseModel()
    {
        try
        {
            // Get host manager reference.
            hostBrowser_ = App.getHostBrowser();

            if (hostBrowser_ != null)
            {
                // Set root node's context info object.
                root_.removeAllChildren();
                root_.setUserObject
                (
                    new ContextInfo("Hosts", ContextInfo.ROOTCONTEXT)
                );

                // Add a host node for each host.
                HostDescriptor[] licensedHosts = hostBrowser_.hosts();
                for (int i=0; i <licensedHosts.length; i++)
                {

                    // Create new contextInfo object for this FoamX host.
                    ContextInfo contextInfo;
                    if (licensedHosts[i].alive)
                    {
                        contextInfo = new ContextInfo
                        (
                            licensedHosts[i].name,
                            ContextInfo.HOST
                        );
                    }
                    else
                    {
                        contextInfo = new ContextInfo
                        (
                            licensedHosts[i].name,
                            ContextInfo.HOSTOFFLINE
                        );
                    }


                    contextInfo.setStatusText
                    (
                        "FoamX Host : " + contextInfo.toString()
                    );

                    // Insert node into tree.
                    DefaultMutableTreeNode node = new DefaultMutableTreeNode
                    (
                        contextInfo,
                        true
                    );

                    insertNodeInto(node, root_, root_.getChildCount());

                    // Add to host map.
                    hostMap_.put(licensedHosts[i].name, node);
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** 
      * Refresh the tree model. Use forceReread only if reread from disk
      */
    public void refresh(boolean forceReread)
    {
        try
        {
            // If we don't yet have a host manager reference,
            // try and re-initialise the model.
            if (hostBrowser_ == null) initialiseModel();

            // Clear map from root/case to treeInfo
            keyMap_.clear();

            // Loop over all host nodes.
            Enumeration iter = hostMap_.keys();
            while (iter.hasMoreElements())
            {
                String hostName = (String)iter.nextElement();

                // Get host node.
                DefaultMutableTreeNode hostNode = 
                    (DefaultMutableTreeNode)hostMap_.get(hostName);
                ContextInfo context = (ContextInfo)hostNode.getUserObject();

                // See if a case browser is running for this host.
                if (context.getCaseBrowser() != null)
                {
                    // Remove all child nodes from the host node.
                    hostNode.removeAllChildren();

                    if (forceReread)
                    {
                        // Refresh caseBrowser's information
                        context.getCaseBrowser().refreshCaseList();
                    }

                    // Re-initialise the hostnode.
                    addCaseBrowserNodes(hostNode);
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** 
      * Refresh the tree model for one root directory only.
      */
    public void refreshRoot(String hostName, String caseRoot)
    {
        try
        {
            // If we don't yet have a host manager reference,
            // try and re-initialise the model.
            if (hostBrowser_ == null) initialiseModel();

            // Find context for give host
            DefaultMutableTreeNode hostNode = 
                (DefaultMutableTreeNode)hostMap_.get(hostName);
            ContextInfo context = (ContextInfo)hostNode.getUserObject();

            // See if a case browser is running for this host.
            if (context.getCaseBrowser() != null)
            {
                if (!rootDirMap_.containsKey(caseRoot))
                {
                    throw new FoamXException
                    (
                        "Invalid caseRoot '" + caseRoot + "'"
                    );
                }

                // Re-initialise the hostnode.
                addRootDirNodes
                (
                    hostNode,
                    caseRoot,
                    (DefaultMutableTreeNode)rootDirMap_.get(caseRoot));
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }


    //--------------------------------------------------------------------------
    /** 
      * Refresh the host tree model.
      */
    public void refreshHosts()
    {
        try
        {
            // If we don't yet have a host manager reference,
            // try and re-initialise the model.
            if (hostBrowser_ == null) initialiseModel();

            // Have hostbrowser determine accessability of hosts
            hostBrowser_.refreshHostList();

            // Update information on tree
            HostDescriptor[] licensedHosts = hostBrowser_.hosts();
            for (int i=0; i <licensedHosts.length; i++)
            {
                // Get host node.
                DefaultMutableTreeNode hostNode = 
                    (DefaultMutableTreeNode)hostMap_.get(licensedHosts[i].name);
                if (hostNode == null)
                {
                    throw new FoamXException
                    (
                        "Host not found in hostMap : " + licensedHosts[i].name
                    );
                }

                ContextInfo context = (ContextInfo)hostNode.getUserObject();

                if (licensedHosts[i].alive)
                {
                    context.setType(ContextInfo.HOST);
                }
                else
                {
                    context.setType(ContextInfo.HOSTOFFLINE);
                }
            }
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (FoamXException fxEx)
        {
            App.handleException(fxEx);
        }
    }

//    //--------------------------------------------------------------------------
//    /** Get hostName for given caseBrowser */
//    public String getHostName(ICaseBrowser caseBrowser)
//    {
//
//        // Loop over all host nodes.
//        itereration iter = hostMap_.keys();
//        while (iter.hasMoreElements())
//        {
//            String hostName = (String)iter.nextElement();
//
//            // Get host node.
//            DefaultMutableTreeNode hostNode = 
//                (DefaultMutableTreeNode)hostMap_.get(hostName);
//            ContextInfo context  = (ContextInfo)hostNode.getUserObject();
//
//            if (context.getCaseBrowser()  == caseBrowser)
//            {
//                return hostName;
//            }
//        }
//        return null;
//    }

    //--------------------------------------------------------------------------

    // Connects to existing/starts new caseBrowser. Returns true if started.
    static private void getCaseBrowserReference
    (
        final String hostName,
        final IHostBrowser hostBrowser,
        ICaseBrowserHolder browserHolder
    ) throws FoamXException, FoamXSYSError, FoamXError, FoamXIOError
    {
        boolean hostOK = false;

        int sleep = Integer.parseInt
        (
            App.getOptions().getProperty("FoamX.Sleep", "500")
        );
        int nRetries = Integer.parseInt
        (
            App.getOptions().getProperty("FoamX.NRetries", "20")
        );

        // First try to connect to already running browser
        hostOK = hostBrowser.getCaseBrowserReference
        (
            hostName,
            browserHolder
        );
        if (hostOK)
        {
            ICaseBrowser caseBrowser = browserHolder.value;
            String[] options = { "Use", "Cancel" };
            int ret =
                JOptionPane.showOptionDialog
                (
                    App.getRootFrame(),
                    "Found running Case Browser"
                  + "\nUse it ?",
                    "Case Browser",
                    JOptionPane.DEFAULT_OPTION,
                    JOptionPane.QUESTION_MESSAGE,
                    null,
                    options,
                    options[0]
                );
            if (ret == 0) return;
        }

        // Start case browser
        hostBrowser.openCaseBrowser(hostName);

        // Go into loop to detect when it has registered.
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

                hostOK = hostBrowser.getCaseBrowserReference
                (
                    hostName,
                    browserHolder
                );

                if (hostOK) return;
            }

            if (!hostOK)
            {
                String[] options = { "Retry", "Exit" };
                int ret =
                    JOptionPane.showOptionDialog
                    (
                        App.getRootFrame(),
                        "Case Browser cannot be contacted.",
                        "Case Browser",
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
                        "No caseBrowser started on "
                      + hostName
                    );
                }
            }
        }
    }


    //--------------------------------------------------------------------------

    boolean startCaseBrowser(String hostName)
    {
        boolean ret = false;

        try
        {
            // Lookup the host node.
            if (!hostMap_.containsKey(hostName))
            {
                throw new FoamXException
                (
                    "Invalid host name '" + hostName + "'"
                );
            }

            DefaultMutableTreeNode hostNode = 
                (DefaultMutableTreeNode)hostMap_.get(hostName);
            ContextInfo context = (ContextInfo)hostNode.getUserObject();

            // Make sure that the case browser is not running on this host.
            if (context.getCaseBrowser() != null)
            {
                throw new FoamXException
                (
                    "Case browser already running on host'" + hostName + "'"
                );
            }

            ICaseBrowserHolder holder = new ICaseBrowserHolder();

            getCaseBrowserReference(hostName, hostBrowser_, holder);

            // Attach the case browser object to the host node.
            ICaseBrowser caseBrowser = holder.value;
            context.setCaseBrowser(caseBrowser);

            // Build the case details under this node.
            addCaseBrowserNodes(hostNode);

            // Return success code.
            ret = true;
        }
        catch (FoamXSYSError sysErr)
        {
            App.handleException(sysErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (FoamXException fxEx)
        {
            App.handleException(fxEx);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return ret;
    }

    //--------------------------------------------------------------------------

    // Change contextInfo to show host is inaccessible
    void setHostOffLine(String hostName, String msg)
    {
        try
        {
            // Lookup the host node.
            if (!hostMap_.containsKey(hostName))
            {
                throw new FoamXException
                (
                    "Invalid host name '" + hostName + "'"
                );
            }

            DefaultMutableTreeNode hostNode = 
                (DefaultMutableTreeNode)hostMap_.get(hostName);
            ContextInfo context = (ContextInfo)hostNode.getUserObject();

            // Mark this server as offline.
            context.setType(ContextInfo.HOSTOFFLINE);
            context.setStatusText(msg);

            // Tell hostBrowser
            hostBrowser_.hostIsDead(hostName);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (FoamXException fxEx)
        {
            App.handleException(fxEx);
        }
    }

    //--------------------------------------------------------------------------

    void stopCaseBrowser(String hostName)
    {
        try
        {
            // Lookup the host node.
            if (!hostMap_.containsKey(hostName))
            {
                throw new FoamXException
                (
                    "Invalid host name '" + hostName + "'"
                );
            }

            DefaultMutableTreeNode hostNode =
                (DefaultMutableTreeNode)hostMap_.get(hostName);

            ContextInfo context = (ContextInfo)hostNode.getUserObject();

            // Make sure that the case browser is running on this host.
            ICaseBrowser caseBrowser = context.getCaseBrowser();
            if (caseBrowser == null)
            {
                throw new FoamXException
                (
                    "Case browser not running on host '" + hostName + "'"
                );
            }

//            // Clear the cases map for cases managed by this caseBrowser.
//            CaseDescriptor[] cases = caseBrowser.cases();
//            for (int i = 0; i <cases.length; i++)
//            {
//                String key =
//                    cases[i].rootDir + "/" + cases[i].caseName;
//
//                if (!keyMap_.containsKey(key))
//                {
//                    throw new FoamXException
//                    (
//                        "Case not in tree : " + key
//                    );
//                }
//                keyMap_.remove(key);
//            }

            keyMap_.clear();

            // Detach the case browser object from the host node.
            context.setCaseBrowser(null);

            // Close the C++ caseBrowser.
            caseBrowser.close();

            // Remove all child nodes from the host node.
            hostNode.removeAllChildren();

            
        }
        catch (FoamXException fxEx)
        {
            App.handleException(fxEx);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    // Build caseRoot nodes for hostNode
    protected void addCaseBrowserNodes(DefaultMutableTreeNode hostNode)
    {
        try
        {
            ContextInfo hostContext = (ContextInfo)hostNode.getUserObject();
            ICaseBrowser caseBrowser = hostContext.getCaseBrowser();

            if (caseBrowser == null)
            {
                throw new FoamXException("Invalid case browser reference.");
            }

            // Get a list of cases served by this case browser.
            String[] roots =
                caseBrowser.foamProperties().rootDirectories();
            String[] rawRoots =
                caseBrowser.foamProperties().rawRootDirectories();

            // Loop over all case roots.
            for (int i = 0; i <roots.length; i++)
            {
                // Add a new case root node.
                // Use original name as display name but expanded one as
                // case root.
                ContextInfo caseRootContextInfo = new ContextInfo
                (
                    caseBrowser,
                    roots[i],
                    rawRoots[i]
                );

                // Add a new case root node.
                caseRootContextInfo.setStatusText
                (
                    "Root Directory : " + rawRoots[i]
                );

                DefaultMutableTreeNode caseRootNode =
                    new DefaultMutableTreeNode
                    (
                        caseRootContextInfo,
                        true
                    );

                insertNodeInto
                (
                    caseRootNode,
                    hostNode,
                    hostNode.getChildCount()
                );

                // Store reference to node
                rootDirMap_.put(roots[i], caseRootNode);
//
//                // Add cases
//                addRootDirNodes(hostNode, roots[i], caseRootNode);
            }
        }
        catch (FoamXException fxEx)
        {
            App.handleException(fxEx);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    // Insert cases below caseRoot node
    protected void addRootDirNodes
    (
        DefaultMutableTreeNode hostNode,
        String caseRoot,
        DefaultMutableTreeNode caseRootNode
    )
    {
        try
        {
            ContextInfo hostContext = (ContextInfo)hostNode.getUserObject();
            ICaseBrowser caseBrowser = hostContext.getCaseBrowser();

            if (caseBrowser == null)
            {
                throw new FoamXException("Invalid case browser reference.");
            }

            // Refresh caseBrowser's information
            caseBrowser.addToCaseList(caseRoot);

            // Loop over all cases and add case nodes.
            CaseDescriptor[] cases = caseBrowser.cases();

            // Clear the caseRoot tree 
            caseRootNode.removeAllChildren();

            for (int i = 0; i < cases.length; i++)
            {
                if (cases[i].rootDir.equals(caseRoot))
                {
                    // Add case node.
                    ContextInfo caseNameContextInfo = new ContextInfo
                    (
                        caseBrowser,
                        cases[i]
                    );

                    caseNameContextInfo.setStatusText
                    (
                        "Case : " + cases[i].caseName
                      + " (" + cases[i].app + ")"
                    );

                    DefaultMutableTreeNode caseNode = new DefaultMutableTreeNode
                    (
                        caseNameContextInfo,
                        false
                    );

                    insertNodeInto
                    (
                        caseNode,
                        caseRootNode,
                        caseRootNode.getChildCount()
                    );

                    // Update key to tree mapping
                    String key = cases[i].rootDir + "/" + cases[i].caseName;
                    keyMap_.put(key, caseNode);
                }
            }
        }
        catch (FoamXException fxEx)
        {
            App.handleException(fxEx);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Updates single case node.
     *  Used for case closing only.
     */
    protected void refreshCaseNode
    (
        String caseRoot,
        String caseName
    )
    {
        try
        {
            // See if we have a treeNode for this case.
            String key = caseRoot + "/" + caseName;

            if (!keyMap_.containsKey(key))
            {
                return;
            }

            DefaultMutableTreeNode caseNode =
                (DefaultMutableTreeNode)keyMap_.get(key);
            ContextInfo caseNameContextInfo =
                (ContextInfo)caseNode.getUserObject();
            ICaseBrowser caseBrowser = caseNameContextInfo.getCaseBrowser();

            // Update contextinfo for caseNode
            // Go back to caseBrowser list and look for caseDescriptor for
            // my case
            //caseBrowser.addToCaseList(caseRoot);
            CaseDescriptor[] cases = caseBrowser.cases();
            CaseDescriptor caseDescriptor = null;
            for (int i = 0; i <cases.length; i++)
            {
                if
                (
                    (caseRoot.equals(cases[i].rootDir))
                 && (caseName.equals(cases[i].caseName))
                )
                {
                    caseDescriptor = cases[i];
                    break;
                }
            }
            if (caseDescriptor == null)
            {
                throw new FoamXException
                (
                    "Case not in caseBrowser : " + key
                );
            }

            // Release old references.
            cases = null;
            caseNameContextInfo.setCaseBrowser(null);
            caseNameContextInfo.setCaseDescriptor(null);

            // Update user object on tree node with current caseDescriptor
            caseNode.setUserObject
            (
                new ContextInfo(caseBrowser, caseDescriptor)
            );
        }
        catch (FoamXException fxEx)
        {
            App.handleException(fxEx);
        }
    }

    //--------------------------------------------------------------------------
    /** Deletes single case node.
     */
    protected void deleteCaseNode
    (
        String caseRoot,
        String caseName
    )
    {
        try
        {
            // See if we have a treeNode for this case.
            String key = caseRoot + "/" + caseName;

            if (!keyMap_.containsKey(key))
            {
                // Not in tree
                throw new FoamXException
                (
                    "Case not in caseBrowser : " + key
                );
            }

            // Remove treeNode
            DefaultMutableTreeNode caseNode =
                (DefaultMutableTreeNode)keyMap_.get(key);
            DefaultMutableTreeNode rootNode =
                (DefaultMutableTreeNode)caseNode.getParent();

            ContextInfo caseNameContextInfo =
                (ContextInfo)caseNode.getUserObject();

            // Release old references.
            caseNameContextInfo.setCaseBrowser(null);
            caseNameContextInfo.setCaseDescriptor(null);

            // Remove caseNode from tree
            rootNode.remove(caseNode);
        }
        catch (FoamXException fxEx)
        {
            App.handleException(fxEx);
        }
    }

    //--------------------------------------------------------------------------
    /** Close all case browsers. */
    public void shutdown()
    {
        // Clear keyMap
        keyMap_.clear();

        // Loop over all host nodes.
        Enumeration iter = hostMap_.keys();
        while (iter.hasMoreElements())
        {
            String hostName = (String)iter.nextElement();

            // Get host node.
            DefaultMutableTreeNode hostNode =
                (DefaultMutableTreeNode)hostMap_.get(hostName);
            ContextInfo context  = (ContextInfo)hostNode.getUserObject();

            // See if a case browser is running on this host.
            ICaseBrowser caseBrowser = context.getCaseBrowser();
            if (caseBrowser != null)
            {
                // Close the case browser on this node.
                caseBrowser.close();
            }
        }
    }
}


//--------------------------------------------------------------------------
