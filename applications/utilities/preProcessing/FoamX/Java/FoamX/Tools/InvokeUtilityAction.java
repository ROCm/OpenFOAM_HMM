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
package FoamX.Tools;

import java.awt.event.ActionEvent;
import java.util.Vector;
import javax.swing.*;

import org.omg.CORBA.StringHolder;

import FoamXServer.IDictionaryEntryHolder;
import FoamXServer.CaseBrowser.ICaseBrowser;
import FoamXServer.CaseServer.ICaseServer;
import FoamXServer.ApplicationDescriptor;

import FoamX.App;
import FoamX.Editors.StringListEditor;
import FoamX.Editors.DictionaryEntryEditor.EntryCache.DictionaryEntryCache;
import FoamX.Editors.DictionaryEntryEditor.EntryCache.DictionaryCache;
import FoamX.Exceptions.FoamXCancel;
import FoamX.Util.BusyCursor;
import FoamX.Util.RunPanel;


public class InvokeUtilityAction
    extends AbstractAction
{
    //--------------------------------------------------------------------------

    protected ICaseBrowser caseBrowser_;
    protected String caseRoot_;
    protected String caseName_;
    protected ICaseServer caseServer_;

    //--------------------------------------------------------------------------
    /** InvokeUtilityAction constructor for a case browser. */
    public InvokeUtilityAction
    (
        ICaseBrowser caseBrowser,
        String caseRoot,
        String caseName
    )
    {
        caseBrowser_ = caseBrowser;
        caseRoot_    = caseRoot;
        caseName_    = caseName;
        caseServer_  = null;

        putValue
        (
            Action.SMALL_ICON,
            App.getResources().getIcon("InvokeUtilityImage")
        );
        putValue(Action.NAME, "Invoke Foam Utility");
        putValue(Action.SHORT_DESCRIPTION, "Invoke Foam Utility");
        putValue(Action.LONG_DESCRIPTION, "Invoke Foam Utility");
    }

//    //--------------------------------------------------------------------------
//    /** InvokeUtilityAction constructor for a case server. */
//    public InvokeUtilityAction(ICaseServer caseServer)
//    {
//        caseBrowser_ = null;
//        caseRoot_    = null;
//        caseName_    = null;
//        caseServer_  = caseServer;
//
//        putValue
//        (
//            Action.SMALL_ICON,
//            App.getResources().getIcon("InvokeUtilityImage")
//        );
//        putValue(Action.NAME, "Invoke Foam Utility");
//        putValue(Action.SHORT_DESCRIPTION, "Invoke Foam Utility");
//        putValue(Action.LONG_DESCRIPTION, "Invoke Foam Utility");
//    }
//
    //--------------------------------------------------------------------------

    public void actionPerformed(ActionEvent evt)
    {
        try
        {
            // Source must have a utilityDescriptor client property.
            if (evt.getSource() instanceof JComponent)
            {
                // Get ApplicationDescriptor object.
                JComponent comp = (JComponent)evt.getSource();
                ApplicationDescriptor desc =
                    (ApplicationDescriptor)comp.getClientProperty
                    (
                        "utilityDescriptor"
                    );

                // Invoke the utility on the specified server.
//                if (caseServer_ != null)
//                {
//                    invokeCaseServerUtility(desc);
//                }
//                else if (caseBrowser_ != null)
//                {
                    invokeCaseBrowserUtility(desc);
//                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    /** construct argument list for utility */
    static public String[] makeUtilityArgs
    (
        ApplicationDescriptor desc,
        String caseRoot,
        String caseName
    )
        throws FoamXCancel
    {
        // Construct argument list.

        String[] caseArgs = null;

//        // If arguments are required, prompt the user for them.
//        if (desc.arguments.length> 0)
//        {
//            Vector labels = new Vector(desc.arguments.length);
//            Vector values = new Vector(desc.arguments.length);
//
//            // Use root, case as first arguments
//
//            if (desc.caseTool)
//            {
//                labels.addElement("Case Root");
//                labels.addElement("Case Name");
//
//                values.addElement(caseRoot);
//                values.addElement(caseName);
//            }
//
//            // Add other arguments as from the toolDescriptor
//
//            for (int i = 0; i <desc.arguments.length; i++)
//            {
//                labels.addElement(desc.arguments[i]);
//                values.addElement(new String());
//            }
//
//            // Pop up panel to edit argument list
//
//            StringListEditor dlg =
//                new StringListEditor
//                (
//                    App.getRootFrame(),
//                    desc.name + " Arguments",
//                    labels, values
//                );
//
//            dlg.show();
//            if (dlg.wasCancelled())
//            {
//                throw new FoamXCancel
//                (
//                    "Editing cancelled"
//                );
//            }
//
//            caseArgs = new String[values.size()];
//
//            for (int i = 0; i <values.size(); i++)
//            {
//                caseArgs[i] = (String)values.elementAt(i);
//            }
//        }
//        else
//        {
//            if (desc.caseTool)
//            {
                caseArgs = new String[2];
                caseArgs[0] = "-case";
                caseArgs[1] = caseRoot + '/' + caseName;
//            }
//            else
//            {
//                caseArgs = null;
//            }
//        }
//
        return caseArgs;
    }

    //--------------------------------------------------------------------------

//    public void invokeCaseServerUtility(ApplicationDescriptor desc)
//    {
//        try
//        {
//            // See if a client side bean has been defined for this action.
//            if (desc.clientBean.length()> 0)
//            {
//                // Load the specified bean.
//                Object utilityObj;
//                try
//                {
//                    ClassLoader loader = ClassLoader.getSystemClassLoader();
//                    utilityObj =
//                        java.beans.Beans.instantiate(loader, desc.clientBean);
//                }
//                catch (Exception ex)
//                {
//                    throw new Exception
//                    (
//                        "Couldn't instantiate utility bean "
//                        + desc.clientBean + "."
//                    );
//                }
//
//                // Make sure that this bean supports the IUtility interface.
//                java.lang.Class c;
//                c = Class.forName("FoamX.Tools.IUtility");
//                if (!java.beans.Beans.isInstanceOf(utilityObj, c))
//                {
//                    throw new Exception
//                    (
//                        "Utility bean " + desc.clientBean
//                        + " does not support the IUtility interface."
//                    );
//                }
//
//                // Invoke the utility.
//                IUtulity utilityBean = (IUtility)utilityObj;
//                utilityBean.prepareUtility();
//                utilityBean.invokeCaseUtility
//                (
//                    desc,
//                    caseServer_,
//                    caseServer_.caseRoot(),
//                    caseServer_.caseName()
//                );
//            }
//            else
//            {
//                // Construct argument list.
//                String[] caseArgs = makeUtilityArgs
//                (
//                    desc,
//                    caseServer_.caseRoot(),
//                    caseServer_.caseName()
//                );
//
//                // Show busy cursor.
//                BusyCursor busyCursor = new BusyCursor(App.getRootFrame());
//
//                // Invoke the utility.
//                try
//                {
//                    caseServer_.invokeUtility(desc.name, caseArgs);
//
//                    JOptionPane.showMessageDialog
//                    (
//                        App.getRootFrame(),
//                        "The Foam utility " + desc.name + " has been executed.",
//                        "Foam Utility",
//                        JOptionPane.PLAIN_MESSAGE
//                    );
//                }
//                catch (Exception ex)
//                {
//                    App.handleAllExceptions(ex);
//                }
//            }
//        }
//        catch (FoamXCancel ex)
//        {
//            // Editing of arguments was cancelled.
//        }
//        catch (Exception ex)
//        {
//            App.handleAllExceptions(ex);
//        }
//    }

    //--------------------------------------------------------------------------

    public void invokeCaseBrowserUtility(ApplicationDescriptor desc)
    {
        try
        {
            // Get the DictionaryEntry reference for this dictionary.
            IDictionaryEntryHolder dictHolder = new IDictionaryEntryHolder();
            caseBrowser_.foamProperties().getUtilityControlDict
            (
                desc.name,
                caseRoot_,
                caseName_,
                dictHolder
            );


            DictionaryCache controlDict = null;
            if (dictHolder.value != null)
            {
                controlDict = (DictionaryCache)DictionaryEntryCache.New
                (
                    dictHolder.value
                );
            }

            RunPanel runPanel = new RunPanel
            (
                App.getRootFrame(),
                caseBrowser_,
                caseRoot_,
                caseName_,
                desc.name,
                false,              // Is utility
                desc.name + "Dict",
                controlDict
            );
            runPanel.setVisible(true);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
