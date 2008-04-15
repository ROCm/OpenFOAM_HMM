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

package FoamX;

import java.awt.*;
import java.awt.event.*;

import java.io.*;
import java.net.URL;
import java.util.*;

import javax.swing.*;

import org.omg.CosNaming.*;

import FoamXServer.FoamXIOError;
import FoamXServer.FoamXError;
import FoamXServer.FoamXSYSError;

import FoamXServer.HostBrowser.IHostBrowserHelper;
import FoamXServer.HostBrowser.IHostBrowser;

import FoamX.Exceptions.FoamXException;
import FoamX.Exceptions.CaseBrowserIOException;
import FoamX.Util.FileEditor;
import FoamX.Util.MemoryView;
import FoamX.Util.MenuItemActionChangedListener;
import FoamX.Util.ResourceHelper;

import FoamX.Options.OptionsManager;
import FoamX.ToolbarManagement.Toolbar;
import FoamX.ToolbarManagement.ToolbarManager;
import FoamX.ActionManagement.ActionManager;
import FoamX.TaskManagement.TaskManager;
import FoamX.Reporting.ReportingManager;
import FoamX.Reporting.ReportingWindow;
import FoamX.CaseManagement.CaseManager;
import FoamX.CaseManagement.CaseBrowserPanel;
import FoamX.WindowManagement.WindowManager;

public class App
    extends javax.swing.JFrame
    implements FoamX.ActionManagement.IActionProvider
{
    //--------------------------------------------------------------------------

    /** Singleton App object */
    protected static App mainFrame_;

    /** Singleton resource helper object. */
    protected static ResourceHelper resources_;

    /** Singleton options object. */
    protected static OptionsManager options_;

    /** Singleton window manager object. */
    protected static WindowManager windowManager_;

    /** Singleton reporting manager object. */
    protected static ReportingManager reportingManager_;

    /** Singleton action manager object. */
    protected static ActionManager actionManager_;

    /** Singleton task manager object. */
    protected static TaskManager taskManager_;

    /** Singleton toolbar manager object. */
    protected static ToolbarManager toolbarManager_;

    /** Singleton Corba host manager object. */
    protected static IHostBrowser hostBrowser_;

    /** Singleton case manager object. */
    protected static CaseManager caseManager_;

    //--------------------------------------------------------------------------

    public static JFrame getRootFrame()
    {
        return mainFrame_;
    }
    public static ResourceHelper getResources()
    {
        return resources_;
    }
    public static OptionsManager getOptions()
    {
        return options_;
    }
    public static WindowManager getWindowManager()
    {
        return windowManager_;
    }
    public static ReportingManager getReportingManager()
    {
        return reportingManager_;
    }
    public static TaskManager getTaskManager()
    {
        return taskManager_;
    }
    public static ActionManager getActionManager()
    {
        return actionManager_;
    }
    public static ToolbarManager getToolbarManager()
    {
        return toolbarManager_;
    }
    public static IHostBrowser getHostBrowser()
    {
        return hostBrowser_;
    }
    public static CaseManager getCaseManager()
    {
        return caseManager_;
    }

    //--------------------------------------------------------------------------

    protected static final String IMAGE_SUFFIX  = "Image";
    protected static final String LABEL_SUFFIX  = "Label";
    protected static final String SHORT_DESCRIPTION_SUFFIX = "ShortDescription";
    protected static final String LONG_DESCRIPTION_SUFFIX = "LongDescription";
    protected static final String ACTION_SUFFIX = "Action";
    protected static final String TIP_SUFFIX    = "Tooltip";

    public static final int DEBUGLEVEL_ERROR   = 0;
    public static final int DEBUGLEVEL_NORMAL  = 1;
    public static final int DEBUGLEVEL_VERBOSE = 2;
    public static final int DEBUGLEVEL_DEBUG   = 3;
    protected static int debugLevel_       = 1;

    protected Toolbar defaultToolbar_;

    /** Actions provided by this object. */
    protected Action[] defaultActions_;

    protected boolean showCallStack_ = true;

    //--------------------------------------------------------------------------
    //---- Static Methods
    //--------------------------------------------------------------------------
    /**
     * Main program entry point.
     * @param args The command line arguments.
     */
    public static void main(String args[])
    {
        mainFrame_ = new App();

        mainFrame_.initialise(args);

        mainFrame_.show();
    }

    //--------------------------------------------------------------------------
    //---- App Methods
    //--------------------------------------------------------------------------
    /** FoamX Application constructor. */
    public App()
    {
        try
        {
            // Initialise the basic GUI.
            initComponents();

            // Sometimes Forte does not add menuBar so do explicitly
            setJMenuBar(menuBar_);


            // Initialise the resource helper object.
            resources_ = new ResourceHelper("resources.FoamX");

            // Initialise the reporting manager object.
            // Any subsequent errors will be reported in the report window.
            reportingManager_ = new ReportingManager();
            horizontalSplitPane_.setBottomComponent(reportingManager_);

            // Initialise the option manager object.
            options_ = new OptionsManager();

            // Initialise the window manager object.
            windowManager_ = new WindowManager(desktopPanel_);

            // Initialise the task status bar.
            taskManager_ = new TaskManager();
            statusPanel_.add(taskManager_.getStatusBar());

            // Create the action manager object.
            actionManager_ = new ActionManager();

            // Create action objects.
            defaultActions_    = new Action[8];
            defaultActions_[0] = new ExitAction("Exit");
            defaultActions_[1] = new ShowMemoryViewAction("ShowMemoryView");
            //defaultActions_[2] = new EditOptionsAction("EditOptions");
            defaultActions_[2] = new ShowHTMLAction("FoamX");
            defaultActions_[3] = new ShowPDFAction("DemoGuide");
            defaultActions_[4] = new ShowPDFAction("ProgrammersGuide");
            defaultActions_[5] = new ShowHTMLAction("Code");
            defaultActions_[6] = new ShowHTMLOnlineAction("Online");
            defaultActions_[7] = new ShowHTMLAction("AboutFoamX");

            // Register all actions with the action manager.
            actionManager_.registerActions(this);
            actionManager_.registerActions(windowManager_);

            // Create toolbar manager object and add to main frame.
            toolbarManager_ = new ToolbarManager();
            GridBagConstraints gbc = new java.awt.GridBagConstraints();
            gbc.fill    = java.awt.GridBagConstraints.HORIZONTAL;
            gbc.gridx   = 0;
            gbc.gridy   = 0;
            gbc.weightx = 1.0;
            gbc.anchor  = java.awt.GridBagConstraints.NORTHWEST;
            getContentPane().add(toolbarManager_, gbc);

            // Initialise the default user interface.
            initialiseDefaultGUI();
            
            //------------------------------------------------------------------
            // Now, initialise the case manager object.
            caseManager_ = new CaseManager();
            verticalSplitPane_.setLeftComponent(caseManager_);

            // Register the reporting manager as a listener to the case browser.
            caseManager_.addCaseManagerListener(reportingManager_);

            //------------------------------------------------------------------
            // Get the debug level.
            debugLevel_ = Integer.parseInt
            (
                options_.getProperty("FoamX.DebugLevel", "1")
            );

            if (debugLevel_ < App.DEBUGLEVEL_ERROR)
            {
                debugLevel_ = App.DEBUGLEVEL_ERROR;
            }

            if (debugLevel_> App.DEBUGLEVEL_DEBUG)
            {
                debugLevel_ = App.DEBUGLEVEL_DEBUG;
            }

            //------------------------------------------------------------------
            // Get the frame size.
            int frameWidth  = Integer.parseInt
            (
                options_.getProperty("FoamX.FrameWidth", "800")
            );

            int frameHeight = Integer.parseInt
            (
                options_.getProperty("FoamX.FrameHeight", "600")
            );

            // Repack and reposition.
            pack();
            java.awt.Dimension screenSize =
                java.awt.Toolkit.getDefaultToolkit().getScreenSize();
            setSize(new java.awt.Dimension(frameWidth, frameHeight));
            setLocation
            (
                Math.max((screenSize.width  - frameWidth  ) / 2, 0),
                Math.max((screenSize.height - frameHeight) / 2, 0)
            );
            horizontalSplitPane_.setDividerLocation(0.90F);
            verticalSplitPane_.setDividerLocation(0.50F);

//            verticalSplitPane_.setSize
//            (
//                new java.awt.Dimension
//                (
//                    300,
//                    verticalSplitPane_.getSize().height
//                )
//            );

//            horizontalSplitPane_.setDividerLocation(0.80F);
//            verticalSplitPane_.setDividerLocation(0.90F);
//            horizontalSplitPane_.setSize
//            (
//                new java.awt.Dimension
//                (
//                    horizontalSplitPane_.getSize().width,
//                    frameHeight - 200
//                )
//            );
        }
        catch (Exception ex)
        {
            // Can't use App's exception handler here
            // since we are in the constructor.
            System.err.println(ex.getMessage());
            ex.printStackTrace();
            System.exit(-1);
        }
    }

    //--------------------------------------------------------------------------
    /** Initialise the FoamX Application object. */
    void initialise(String[] args)
    {
        try
        {
            // Initialise the host manager reference.
            if (initialiseHostBrowser(args))
            {
                // Refresh the case browser to show the hosts.
                caseManager_.getCaseBrowser().refreshCaseBrowser(false);

                // Do some simple command line parsing

                int argI = 0;

                while
                (
                    (argI < args.length)
                 && (args[argI].substring(0,4).equals("-ORB"))
                )
                {
                    argI += 2;
                }

                String hostName = null;
                String caseRoot = null;
                String caseName = null;
                String action = "open";

                if (argI < args.length)
                {
                     hostName = args[argI++];
                }
                if (argI < args.length)
                {
                     caseRoot = args[argI++];
                }
                if (argI < args.length)
                {
                     caseName = args[argI++];
                }
                if (argI < args.length)
                {
                     action = args[argI++];
                }

                if (hostName != null)
                {
                    if
                    (
                        !caseManager_.getCaseBrowser().serverStart
                        (
                            hostName, caseRoot, caseName, action
                        )
                    )
                    {
                        throw new FoamXException
                        (
                            "Cannot find host " + hostName + " in name service"
                        );
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

    public boolean initialiseHostBrowser(String[] args)
    {
        boolean bret = false;

        try
        {
            if (hostBrowser_ == null)
            {
                // Create th'orb.
                org.omg.CORBA.ORB orb = org.omg.CORBA.ORB.init
                (
                    args,
                    options_
                );

                // Get the root naming context from the naming service.
                NamingContextExt     rootContext = null;
                org.omg.CORBA.Object nsObj       = null;

                try
                {
                    // Try resolve_initial_references first.
                    nsObj = orb.resolve_initial_references("PNameService");

                    // Print nice status message.
                    App.printMessage
                    (
                        App.DEBUGLEVEL_VERBOSE,
                        "Connecting to name service "
                      + orb.toString()
                    );
                }
                catch (Exception ex)
                {
                    // The resolve_initial_references call didn't work so try
                    // and read the name service IOR written to the ns.ref file
                    // in the FoamX user directory.
                    String userPath = options_.foamXUserConfigPath();
                    if (userPath != null)
                    {
                        String nsFileName = userPath + "/ns.ref";

                        // Print nice status message.
                        App.printMessage
                        (
                            App.DEBUGLEVEL_VERBOSE,
                            "Reading name service reference from '"
                          + nsFileName + "'."
                        );

                        // Check that the nameservice reference file exists.
                        File nsFile = new File(nsFileName);
                        if (nsFile.exists())
                        {
                            BufferedReader br = new BufferedReader
                            (
                                new FileReader(nsFileName)
                            );
                            String ior = br.readLine();
                            nsObj = orb.string_to_object(ior);

                            // Print nice status message.
                            App.printMessage
                            (
                                App.DEBUGLEVEL_DEBUG,
                                "Connecting to name service."
                            );
                        }
                        else
                        {
                            // Print error message.
                            App.printMessage
                            (
                                App.DEBUGLEVEL_ERROR,
                                "Name service reference file not found in '"
                              + nsFileName + "'."
                            );
                        }
                    }
                }

                try
                {
                    rootContext = NamingContextExtHelper.narrow(nsObj);

                    if (rootContext == null || rootContext._non_existent())
                    {
                        rootContext = null;
                        throw new FoamXException
                        (
                            "Cannot find rootContext in name service."
                        );
                    }

                    // Print success message.
                    App.printMessage
                    (
                        App.DEBUGLEVEL_NORMAL,
                        "Successfully connected to name service."
                    );
                }
                catch (Exception ex)
                {
                    throw new FoamXException
                    (
                        "Default name service could not be contacted."
                    );
                }


                // Get a reference to the host browser.
                try
                {
                    NameComponent[] name = new NameComponent[1];
                    name[0] = new NameComponent("FoamXHostBrowser", "");

                    org.omg.CORBA.Object obj = rootContext.resolve(name);
                    hostBrowser_ = IHostBrowserHelper.narrow(obj);
                    if 
                    (
                        hostBrowser_ == null
                     || hostBrowser_._non_existent()
                    )
                    {
                        hostBrowser_ = null;
                        throw new FoamXException
                        (
                            "Cannot find FoamXHostBrowser in name service."
                        );
                    }
                }
                catch (Exception ex)
                {
                    throw new FoamXException
                    (
                        "Foam Host Browser object could not be contacted."
                    );
                }
            }

            // Indicate success.
            bret = true;
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
            bret = false;
        }

        return bret;
    }

    //--------------------------------------------------------------------------

    public static void handleException(FoamXIOError ioErr)
    {
        // Beep.
        Toolkit.getDefaultToolkit().beep();

        String message =
            ioErr.errorMessage
          + "\nin file " + ioErr.ioFileName + "\n"
          + " start at line " + new Integer(ioErr.ioStartLineNo).toString();
       
        if (ioErr.ioEndLineNo > 0)
        {
           message += " ending at line " + new Integer(ioErr.ioEndLineNo).toString();
        }

        if (debugLevel_ >= App.DEBUGLEVEL_DEBUG)
        {
            message +=
                "\nin function " + ioErr.methodName
              + "\nin file " + ioErr.fileName
              + " at line " + new Integer(ioErr.lineNo).toString();
        }

        // Print the exception's message
        printMessage(message);


        if (debugLevel_ >= App.DEBUGLEVEL_DEBUG)
        {
            ioErr.printStackTrace();
        }

        String fileName = ioErr.ioFileName;
        int subPos = fileName.indexOf("::");
        if (subPos != -1)
        {
            fileName = fileName.substring(0, subPos);
        }
        message += "\n\nDo you want to edit " + fileName + "?";

        // Ask whether to open editor on file
        if
        (
            JOptionPane.showConfirmDialog
            (
                App.getRootFrame(),
                message,
                "FoamXIOError",
                JOptionPane.YES_NO_OPTION,
                JOptionPane.ERROR_MESSAGE
            )
            == JOptionPane.OK_OPTION
        )
        {
            new FileEditor
            (
                true,                   // modal
                false,                  // no follow mode
                null,                   // no caseBrowser
                fileName,
                ioErr.ioStartLineNo,
                ioErr.ioEndLineNo,
                0
            );
        }
    }

    //--------------------------------------------------------------------------

    /* exception same as FoamXIOError but carries extra caseBrowser field */
    public static void handleException(CaseBrowserIOException cbErr)
    {
        FoamXIOError ioErr = cbErr.getIoErr();

        // Beep.
        Toolkit.getDefaultToolkit().beep();

        String message =
            ioErr.errorMessage
          + "\nin file " + ioErr.ioFileName + "\n"
          + " start at line " + new Integer(ioErr.ioStartLineNo).toString();
       
        if (ioErr.ioEndLineNo > 0)
        {
           message += " ending at line " + new Integer(ioErr.ioEndLineNo).toString();
        }

        if (debugLevel_ >= App.DEBUGLEVEL_DEBUG)
        {
            message +=
                "\nin function " + ioErr.methodName
              + "\nin file " + ioErr.fileName
              + " at line " + new Integer(ioErr.lineNo).toString();
        }

        // Print the exception's message
        printMessage(message);

        if (debugLevel_ >= App.DEBUGLEVEL_DEBUG)
        {
            cbErr.printStackTrace();
        }

        String fileName = ioErr.ioFileName;
        int subPos = fileName.indexOf("::");
        if (subPos != -1)
        {
            fileName = fileName.substring(0, subPos);
        }
        message += "\n\nDo you want to edit " + fileName + "?";

        // Ask whether to open editor on file
        if
        (
            JOptionPane.showConfirmDialog
            (
                App.getRootFrame(),
                message,
                "FoamXIOError",
                JOptionPane.YES_NO_OPTION,
                JOptionPane.ERROR_MESSAGE
            )
            == JOptionPane.OK_OPTION
        )
        {
            new FileEditor
            (
                true,                   // modal
                false,                  // edit, not follow mode
                cbErr.getCaseBrower(),  // use caseBrowser to get file
                fileName,
                ioErr.ioStartLineNo,
                ioErr.ioEndLineNo,
                0
            );
        }
    }

    //--------------------------------------------------------------------------

    public static void handleException(FoamXError fxErr)
    {
        // Beep.
        Toolkit.getDefaultToolkit().beep();

        String message = fxErr.errorMessage;

        if (debugLevel_ >= App.DEBUGLEVEL_DEBUG)
        {
            message +=
                "\nin function " + fxErr.methodName
              + "\nin file " + fxErr.fileName
              + " at line " + new Integer(fxErr.lineNo).toString();
        }

        // Print the exception's message
        printMessage(message);

        //TODO: originates maybe from opening case so reportingWindow not valid
        //printStack(fxErr);
        if (debugLevel_ >= App.DEBUGLEVEL_DEBUG)
        {
            fxErr.printStackTrace();
        }

        // Print the message in a dialogue box
        JOptionPane.showMessageDialog
        (
            App.getRootFrame(),
            message,
            "FoamXError",
            JOptionPane.ERROR_MESSAGE
        );
    }

    //--------------------------------------------------------------------------

    /** machine can't be reached */
    public static void handleException(FoamXSYSError sysErr)
    {
        // Beep.
        Toolkit.getDefaultToolkit().beep();

        String message = sysErr.errorMessage;

        if (debugLevel_ >= App.DEBUGLEVEL_DEBUG)
        {
            message +=
                "\nin function " + sysErr.methodName
              + "\nin file " + sysErr.fileName
              + " at line " + new Integer(sysErr.lineNo).toString();
        }

        // Print the exception's message
        printMessage(message);

        // Print the message in a dialogue box
        JOptionPane.showMessageDialog
        (
            App.getRootFrame(),
            message,
            "FoamXSYSError",
            JOptionPane.ERROR_MESSAGE
        );

        // Refresh the case browser to show the hosts.
        CaseBrowserPanel caseBrowserPanel = caseManager_.getCaseBrowser();
        if (caseBrowserPanel != null)
        {
            caseBrowserPanel.setHostOffLine
            (
                sysErr.hostName,
                ""
            );
        }
        //caseBrowserPanel.refreshCaseBrowser(false);
    }

    //--------------------------------------------------------------------------

    public static void handleException(FoamXException fxEx)
    {
        // Beep.
        Toolkit.getDefaultToolkit().beep();

        //printStack(fxErr);
        fxEx.printStackTrace();

        // Print the exception's message
        printMessage(fxEx.getMessage());

        // Print the message in a dialogue box
        JOptionPane.showMessageDialog
        (
            App.getRootFrame(),
            fxEx.getMessage(),
            "FoamXException",
            JOptionPane.ERROR_MESSAGE
        );
    }

    //--------------------------------------------------------------------------

    public static void handleException(MissingResourceException ex)
    {
        // Beep.
        Toolkit.getDefaultToolkit().beep();

        //printStack(fxErr);
        ex.printStackTrace();


        // Print the exception's message
        printMessage(ex.getMessage());

        // Print the message in a dialogue box
        JOptionPane.showMessageDialog
        (
            App.getRootFrame(),
            ex.getMessage()
          + " for key " + ex.getKey() + " for resource " + ex.getClassName(),
            "MissingResourceException",
            JOptionPane.ERROR_MESSAGE
        );
    }

    //--------------------------------------------------------------------------

    public static void handleException(Exception ex)
    {
        // Beep.
        Toolkit.getDefaultToolkit().beep();

        //printStack(fxErr);
        ex.printStackTrace();

        // Print the exception's message
        printMessage(ex.getMessage());

        // Print the message in a dialogue box
        JOptionPane.showMessageDialog
        (
            App.getRootFrame(),
            ex.getMessage(),
            ex.getClass().toString(),
            JOptionPane.ERROR_MESSAGE
        );
    }

    //--------------------------------------------------------------------------

    public static void handleAllExceptions(Exception ex)
    {
        if (ex instanceof FoamXIOError)
        {
            handleException((FoamXIOError)ex);
        }
        else if (ex instanceof FoamXError)
        {
            handleException((FoamXError)ex);
        }
        else if (ex instanceof FoamXSYSError)
        {
            handleException((FoamXSYSError)ex);
        }
        else if (ex instanceof FoamXException)
        {
            handleException((FoamXException)ex);
        }
        else if (ex instanceof CaseBrowserIOException)
        {
            handleException((CaseBrowserIOException)ex);
        }
        else if (ex instanceof MissingResourceException)
        {
            handleException((MissingResourceException)ex);
        }
        else
        {
            handleException(ex);
        }
    }

    //--------------------------------------------------------------------------

    public static void printMessage(String message)
    {
        // Print message to reporting window if we have one.
        if (reportingManager_ != null)
        {
            ReportingWindow reportWin = reportingManager_.getReportingWindow
            (
                "FoamX"
            );

            if (reportWin != null)
            {
                // Print message into main reporting window.
                reportWin.printMessage(message);
            }
            else
            {
                System.err.println(message);
            }
        }
        else
        {
            System.err.println(message);
        }
    }

    //--------------------------------------------------------------------------

    public static void printStack(Exception ex)
    {
        // Print message to reporting window if we have one.
        if (reportingManager_ != null)
        {
            ReportingWindow reportWin = reportingManager_.getReportingWindow
            (
                "FoamX"
            );

            if (reportWin != null)
            {
                ByteArrayOutputStream strStream = new ByteArrayOutputStream();
                PrintStream printStream = new PrintStream(strStream);

                // Print message into main reporting window.
                ex.printStackTrace(printStream);
            }
            else
            {
                // Print the stack trace to standard error.
                ex.printStackTrace();
            }
        }
    }

    //--------------------------------------------------------------------------

    public static void printMessage(int level, String message)
    {
        if (debugLevel_ >= level)
        {
            printMessage(message);
        }
    }


    //--------------------------------------------------------------------------

    public static void pm(Object o)
    {
	java.lang.reflect.Method ms[] = o.getClass().getMethods();
	System.out.println(o.getClass().getName());
	for(int i = 0; i < ms.length; i++)
        {
            System.out.print
            (
                ms[i].getReturnType().getName()
              + " "
              + ms[i].getName()
              + "("
            );
            java.lang.Class cs[] = ms[i].getParameterTypes();
            for(int j = 0; j < cs.length; j++)
            {
                System.out.println
                (
                    "  " + cs[j].getName() +
                    ((j==cs.length-1)?"":",")
                );
            }
            System.out.println(")");
        }
    }

    //--------------------------------------------------------------------------

    protected void initialiseDefaultGUI()
    {
        try
        {
            // Initialise menu and tool bars.
            initialiseMenubar();
            initialiseToolbar();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Initialise the menubar for the app. By default this pulls the
     * definition of the menu from the associated resource file.
     */
    protected void initialiseMenubar()
    {
        try
        {
            menuBar_.removeAll();
            
            // Get names of menu items to create from resource file.
            String[] menuKeys =
                resources_.getResourceStringArray("FoamX.Menubar");

            for (int i=0; i <menuKeys.length; i++)
            {
                // Create a new menu and initialise it.
                JMenu menu = createMenu(menuKeys[i]);
                if (menu != null)
                {
                    menuBar_.add(menu);
                }
            }

            // Add application defined menus.

        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Create a new menu. The menu is defined within
     * the associated resource file.
     */
    protected JMenu createMenu(String key)
    {
        JMenu menu = null;

        try
        {
            String[] itemKeys  = resources_.getResourceStringArray(key);
            String menuLabel =resources_.getResourceString(key + LABEL_SUFFIX);

            // Create menu.
            menu = new JMenu(menuLabel);

            // Use same font as the menu bar.
            menu.setFont(menuBar_.getFont());

            // Add menu items.
            for (int i=0; i <itemKeys.length; i++)
            {
                if (itemKeys[i].equals("-"))
                {
                    // Add separator.
                    menu.addSeparator();
                }
                else
                {
                    // Add menu item.
                    JMenuItem menuItem = createMenuItem(itemKeys[i]);
                    if (menuItem != null)
                    {
                        menu.add(menuItem);
                    }
                }
            }
            menu.setFont(menuBar_.getFont());
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return menu;
    }

    //--------------------------------------------------------------------------
    /**
     * This is the hook through which all menu items are
     * created. It registers the result with the menuitem
     * hashtable so that it can be fetched with getMenuItem().
     */
    protected JMenuItem createMenuItem(String key)
    {
        JMenuItem menuItem = null;

        try
        {
            // Get resource info.
            String menuLabel = resources_.getResourceString
            (
                key + LABEL_SUFFIX
            );
            String actionName = resources_.getResourceString
            (
                key + ACTION_SUFFIX
            );
            URL url = resources_.getResourceURL(key + IMAGE_SUFFIX);

            // Set the action name.
            if (actionName == null) actionName = key;

            // Check that it has been registered.
            if (!actionManager_.isRegistered(actionName))
            {
                throw new Exception
                (
                    "Invalid action for menu item " + key
                  + " in createMenuItem method call.");
            }
            Action action = actionManager_.getAction(actionName);

            // Create and initialise the menu item.
            menuItem = new JMenuItem(menuLabel);

            // Set the icon.
            if (url != null)
            {
                menuItem.setHorizontalTextPosition(JButton.RIGHT);
                menuItem.setIcon(new ImageIcon(url));
            }

            // Attach the specified action.
            if (action != null)
            {
                // Attach action and link up change listener object.
                menuItem.setActionCommand(actionName);
                menuItem.addActionListener(action);
                action.addPropertyChangeListener
                (
                    new MenuItemActionChangedListener(menuItem)
                );
                menuItem.setEnabled(action.isEnabled());
            }
            else
            {
                // Action not found - disable this item.
                menuItem.setEnabled(false);
            }

            // Use same font as the menu bar.
            menuItem.setFont(menuBar_.getFont());
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return menuItem;
    }

    //--------------------------------------------------------------------------
    /**
     * Create the toolbar. By default this reads the
     * resource file for the definition of the toolbar.
     */
    protected void initialiseToolbar()
    {
        try
        {
            // Add default toolbar buttons.
            String[] toolKeys =
                resources_.getResourceStringArray("FoamX.Toolbar");

            if (toolKeys.length > 0)
            {

                // Create default toolbar.
                defaultToolbar_ = toolbarManager_.addToolbar
                (
                    "Default",
                    "Default"
                );


                for (int i=0; i <toolKeys.length; i++)
                {
                    if (toolKeys[i].equals("-"))
                    {
                        defaultToolbar_.addSeparator();
                    }
                    else
                    {
                        addToolbarButton(toolKeys[i]);
                    }
                }
                defaultToolbar_.add(Box.createHorizontalGlue());
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Create a button to go inside of the toolbar. By default this
     * will load an image resource. The image filename is relative to
     * the classpath (including the '.' directory if its a part of the
     * classpath), and may either be in a JAR file or a separate file.
     *
     * @param key The key in the resource file to serve as the basis
     *  of lookups.
     */
    protected void addToolbarButton(String key)
    {
        try
        {
            String  actionName =
                resources_.getResourceString(key + ACTION_SUFFIX);
            URL url = resources_.getResourceURL(key + IMAGE_SUFFIX);
            String tipText = resources_.getResourceString(key + TIP_SUFFIX);
            JButton button = null;

            if (url != null)
            {
                button = new JButton(new ImageIcon(url));

                /*
                button.setRequestFocusEnabled(false);
                button.setMargin(new Insets(1, 1, 1, 1));
                button.setBorder
                (
                    new javax.swing.border.EmptyBorder(1, 1, 1, 1)
                );
                */
                if (actionName == null) actionName = key;

                if (actionManager_.isRegistered(actionName))
                {
                    Action action = actionManager_.getAction(actionName);
                    button.setActionCommand(actionName);
                    button.addActionListener(action);
                }
                else
                {
                    button.setEnabled(false);
                }

                if (tipText != null)
                {
                    button.setToolTipText(tipText);
                }

                defaultToolbar_.add(button);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to initialize the
     *  form.  WARNING: Do NOT modify this code. The content of this method is*
     *  always regenerated by the FormEditor.
     */
    private void initComponents()//GEN-BEGIN:initComponents
    {
        menuBar_ = new javax.swing.JMenuBar();
        jMenu2 = new javax.swing.JMenu();
        jMenuItem2 = new javax.swing.JMenuItem();
        contextMenu_ = new javax.swing.JPopupMenu();
        mnuProperties = new javax.swing.JMenuItem();
        horizontalSplitPane_ = new javax.swing.JSplitPane();
        verticalSplitPane_ = new javax.swing.JSplitPane();
        desktopPanel_ = new javax.swing.JDesktopPane();
        statusPanel_ = new javax.swing.JPanel();
        
        jMenu2.setText("Menu");
        jMenuItem2.setFont(new java.awt.Font("Courier", 0, 10));
        jMenuItem2.setText("Item");
        jMenu2.add(jMenuItem2);
        menuBar_.add(jMenu2);
        contextMenu_.setFont(new java.awt.Font("Courier", 0, 10));
        mnuProperties.setFont(new java.awt.Font("Courier", 0, 10));
        mnuProperties.setText("Properties...");
        contextMenu_.add(mnuProperties);
        
        getContentPane().setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        setTitle("FoamX");
        setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
        setName("MainFrame");
        setBackground(java.awt.Color.white);
        setFont(new java.awt.Font("Courier", 0, 10));
        addWindowListener(new java.awt.event.WindowAdapter()
        {
            public void windowClosing(java.awt.event.WindowEvent evt)
            {
                exitForm(evt);
            }
            public void windowActivated(java.awt.event.WindowEvent evt)
            {
                OnActivate(evt);
            }
        });
        
        horizontalSplitPane_.setOrientation(javax.swing.JSplitPane.VERTICAL_SPLIT);
        horizontalSplitPane_.setOneTouchExpandable(true);
        verticalSplitPane_.setLastDividerLocation(-2);
        verticalSplitPane_.setOneTouchExpandable(true);
        desktopPanel_.setPreferredSize(new java.awt.Dimension(300, 400));
        verticalSplitPane_.setRightComponent(desktopPanel_);
        
        horizontalSplitPane_.setLeftComponent(verticalSplitPane_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        getContentPane().add(horizontalSplitPane_, gridBagConstraints1);
        
        statusPanel_.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 5, 2));
        
        statusPanel_.setBorder(new javax.swing.border.EtchedBorder());
        statusPanel_.setMinimumSize(new java.awt.Dimension(10, 25));
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 2;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.WEST;
        gridBagConstraints1.weightx = 1.0;
        getContentPane().add(statusPanel_, gridBagConstraints1);
        
        pack();
        java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        setSize(new java.awt.Dimension(800, 600));
        setLocation((screenSize.width-800)/2,(screenSize.height-600)/2);
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

  private void OnActivate(java.awt.event.WindowEvent evt)
    {//GEN-FIRST:event_OnActivate
        // Add your handling code here:




  }//GEN-LAST:event_OnActivate

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    /** Exit the Application */
    private void exitForm(java.awt.event.WindowEvent evt)
    {//GEN-FIRST:event_exitForm

        // Invoke the exit action.
        ExitAction exitAction = new ExitAction("Exit");
        exitAction.actionPerformed(null);

    }//GEN-LAST:event_exitForm

    //--------------------------------------------------------------------------

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JMenuBar menuBar_;
    private javax.swing.JMenu jMenu2;
    private javax.swing.JMenuItem jMenuItem2;
    private javax.swing.JPopupMenu contextMenu_;
    private javax.swing.JMenuItem mnuProperties;
    private javax.swing.JSplitPane horizontalSplitPane_;
    private javax.swing.JSplitPane verticalSplitPane_;
    private javax.swing.JDesktopPane desktopPanel_;
    private javax.swing.JPanel statusPanel_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    // ---- IActionProvider interface functions
    //--------------------------------------------------------------------------

    public String getPrefix()
    {
        return "FoamX";
    }

    public javax.swing.Action[] getActions()
    {
        return defaultActions_;
    }

    //--------------------------------------------------------------------------
    //---- Action Classes
    //--------------------------------------------------------------------------

    class ExitAction
        extends AbstractAction
    {
        ExitAction(String ActionResourceName)
        {
            putValue(Action.NAME, ActionResourceName);
            putValue
            (
                Action.SHORT_DESCRIPTION,
                options_.getProperty
                (
                    ActionResourceName + SHORT_DESCRIPTION_SUFFIX, ""
                )
            );
            putValue
            (
                Action.LONG_DESCRIPTION,
                options_.getProperty
                (
                    ActionResourceName + LONG_DESCRIPTION_SUFFIX, ""
                )
            );
        }

        public void actionPerformed(ActionEvent evt)
        {
            if
            (
                JOptionPane.showConfirmDialog
                (
                    App.getRootFrame(),
                    "Are you sure you want to quit FoamX?",
                    "Confirm Quit Operation",
                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.QUESTION_MESSAGE
                ) == JOptionPane.OK_OPTION
            )
            {
                // Unregister the reporting manager as a listener to the case
                // browser object and shut it down.
                caseManager_.removeCaseManagerListener(reportingManager_);
                caseManager_.shutdown();

                // Shutdown the reporting manager object.
                reportingManager_.shutdown();

                // Exit the application.
                class ExitEvent
                    implements Runnable
                {
                    public void run()
                    {
                        System.exit(0);
                    }
                }

                // Put an exit application event into the event dispatch
                // thread event queue.
                SwingUtilities.invokeLater(new ExitEvent());

                //System.exit(0);
            }
        }
    }

    //--------------------------------------------------------------------------

    class ShowMemoryViewAction
        extends AbstractAction
    {
        ShowMemoryViewAction(String ActionResourceName)
        {
            putValue(Action.NAME, ActionResourceName);
            putValue
            (
                Action.SHORT_DESCRIPTION,
                options_.getProperty
                (
                    ActionResourceName + SHORT_DESCRIPTION_SUFFIX, ""
                )
            );
            putValue
            (
                Action.LONG_DESCRIPTION,
                options_.getProperty
                (
                    ActionResourceName + LONG_DESCRIPTION_SUFFIX, ""
                )
            );
        }

        public void actionPerformed(ActionEvent evt)
        {
            MemoryView mv = new MemoryView();
            Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
            Dimension m = mv.getSize ();
            d.width  -= m.width;
            d.height -= m.height;
            d.width  /= 2;
            d.height /= 2;
            mv.setLocation(d.width, d.height);
            mv.setVisible(true);
        }
    }

    //--------------------------------------------------------------------------

//    class EditOptionsAction
//        extends AbstractAction
//    {
//        EditOptionsAction(String ActionResourceName)
//        {
//            putValue(Action.NAME, ActionResourceName);
//            putValue
//            (
//                Action.SHORT_DESCRIPTION,
//                options_.getProperty
//                (
//                    ActionResourceName + SHORT_DESCRIPTION_SUFFIX, ""
//                )
//            );
//            putValue
//            (
//                Action.LONG_DESCRIPTION,
//                options_.getProperty
//                (
//                    ActionResourceName + LONG_DESCRIPTION_SUFFIX, ""
//                )
//            );
//        }
//
//        public void actionPerformed(ActionEvent evt)
//        {
//            options_.showEditDialog(App.getRootFrame());
//        }
//    }

    //--------------------------------------------------------------------------

    class ShowHTMLAction
        extends AbstractAction
    {
        private String HTMLFileResourceName_;

        ShowHTMLAction(String HTMLFileResourceName)
        {
            HTMLFileResourceName_ = HTMLFileResourceName;

            putValue(Action.NAME, "Show" + HTMLFileResourceName_ + "HTML");
            putValue
            (
                Action.SHORT_DESCRIPTION,
                options_.getProperty
                (
                    HTMLFileResourceName_ + SHORT_DESCRIPTION_SUFFIX, ""
                )
            );
            putValue
            (
                Action.LONG_DESCRIPTION,
                options_.getProperty
                (
                    HTMLFileResourceName_ + LONG_DESCRIPTION_SUFFIX, ""
                )
            );
        }

        public void actionPerformed(ActionEvent evt)
        {
            try
            {
                // Construct the default Url of the documentation home page.
                String docRoot = 
                    "file:" + options_.foamXSystemPath() + "/doc";

                // If the URL is specified in FoamXClient.cfg file use that
                docRoot = options_.getProperty("FoamX.DocRoot", docRoot);

                String HTMLFileName = resources_.getResourceString
                (
                    HTMLFileResourceName_ + "HTML"
                );

                // Point the browser object to the HTMLFileName page.
                URL url = new URL(docRoot + "/" + HTMLFileName);

                // Get name of browser.
                String browser = options_.getProperty
                (
                    "FoamX.Browser",
                    "mozilla"
                );

                // Execute the desired browser passing it the required URL.
                Runtime.getRuntime().exec(browser + " " + url.toExternalForm());
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    class ShowHTMLOnlineAction
        extends AbstractAction
    {
        private String HTMLFileResourceName_;

        ShowHTMLOnlineAction(String HTMLFileResourceName)
        {
            HTMLFileResourceName_ = HTMLFileResourceName;

            putValue(Action.NAME, "Show" + HTMLFileResourceName_ + "HTML");
            putValue
            (
                Action.SHORT_DESCRIPTION,
                options_.getProperty
                (
                    HTMLFileResourceName_ + SHORT_DESCRIPTION_SUFFIX, ""
                )
            );
            putValue
            (
                Action.LONG_DESCRIPTION,
                options_.getProperty
                (
                    HTMLFileResourceName_ + LONG_DESCRIPTION_SUFFIX, ""
                )
            );
        }

        public void actionPerformed(ActionEvent evt)
        {
            try
            {
                String HTMLFileName = resources_.getResourceString
                (
                    HTMLFileResourceName_ + "HTML"
                );

                // Point the browser object to the HTMLFileName page.
                URL url = new URL(HTMLFileName);

                // Get name of browser.
                String browser = options_.getProperty
                (
                    "FoamX.Browser",
                    "mozilla"
                );

                // Execute the desired browser passing it the required URL.
                Runtime.getRuntime().exec(browser + " " + url.toExternalForm());
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    class ShowPDFAction
        extends AbstractAction
    {
        private String PDFFileResourceName_;

        ShowPDFAction(String PDFFileResourceName)
        {
            PDFFileResourceName_ = PDFFileResourceName;

            putValue(Action.NAME, "Show" + PDFFileResourceName_ + "PDF");
            putValue
            (
                Action.SHORT_DESCRIPTION,
                options_.getProperty
                (
                    PDFFileResourceName_ + SHORT_DESCRIPTION_SUFFIX, ""
                )
            );
            putValue
            (
                Action.LONG_DESCRIPTION,
                options_.getProperty
                (
                    PDFFileResourceName_ + LONG_DESCRIPTION_SUFFIX, ""
                )
            );
        }

        public void actionPerformed(ActionEvent evt)
        {
            try
            {
                // Construct the default Url of the documentation home page.
                String docRoot = options_.foamXSystemPath() + "/doc";

                // If the URL is specified in FoamXClient.cfg file use that
                docRoot = options_.getProperty("FoamX.DocRoot", docRoot);

                String PDFFileName = resources_.getResourceString
                (
                    PDFFileResourceName_ + "PDF"
                );

                // Get name of browser.
                String browser = options_.getProperty
                (
                    "FoamX.PDFBrowser",
                    "acroread"
                );

System.out.println("Trying to start " + browser + " " + docRoot + "/" + PDFFileName);

                // Execute the desired browser passing it the required file name
                Runtime.getRuntime().exec
                (
                    browser + " " + docRoot + "/" + PDFFileName
                );
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }
}


//------------------------------------------------------------------------------
