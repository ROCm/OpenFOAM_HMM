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

package FoamX.Modules.RawDictionaryEditor;

import java.beans.*;
import java.util.*;
import javax.swing.*;
import java.awt.event.*;

import FoamXServer.*;
import FoamXServer.CaseServer.*;

import FoamX.App;
import FoamX.Util.FoamXTypeEx;
import FoamX.CaseManagement.IModule;
import FoamX.CaseManagement.IModuleHost;
import FoamX.WindowManagement.FoamXInternalFrame;

public class RawDictionaryEditorModule
    extends Object
    implements IModule
{
    //--------------------------------------------------------------------------

    private static final String moduleName_ = "RawDictionaryEditor";
    private static final String moduleDescription_ =
        "Edits Dictionaries in the raw";

    private IModuleHost moduleHost_;
    private int moduleID_;
    private boolean initialised_;
    // Map of dictionary descriptor objects.
    private java.util.Hashtable dictionaryMap_;
    private PropertyChangeSupport propertySupport_;

    //--------------------------------------------------------------------------
    // Module actions.

    public static final String OPENRAWDICTIONARY_ACTION  = "OpenRawDictionary";
    public static final String CLOSERAWDICTIONARY_ACTION = "CloseRawDictionary";
    public static final String SAVERAWDICTIONARY_ACTION  = "SaveRawDictionary";

    public static final int OPENRAWDICTIONARY_INDEX  = 0;
    public static final int CLOSERAWDICTIONARY_INDEX = 1;
    public static final int SAVERAWDICTIONARY_INDEX = 2;
    private Action[] moduleActions_ =
        {
            new OpenRawDictionaryAction(),
            new CloseRawDictionaryAction(),
            new SaveRawDictionaryAction()
        };

    //--------------------------------------------------------------------------
    // Inner class for storing dictionary information.

    public class DictionaryDescriptor
    {
        private String dictionaryName_;
        private String dictionaryPath_;
        private DictionaryTextEditor editor_;
        private int windowID_;

        public DictionaryDescriptor
        (
            String dictionaryName,
            String dictionaryPath
        )
        {
            dictionaryName_ = dictionaryName;
            dictionaryPath_ = dictionaryPath;
            editor_         = null;
            windowID_       = -1;
        }

        public String getName()
        {
            return dictionaryName_;
        }
        public String getPath()
        {
            return dictionaryPath_;
        }
        public DictionaryTextEditor getTextEditor()
        {
            return editor_;
        }
        public void setTextEditor(DictionaryTextEditor editor)
        {
            editor_ = editor;
        }
        public int getWindowID()
        {
            return windowID_;
        }
        public void setWindowID(int windowID)
        {
            windowID_ = windowID;
        }
    }

    //--------------------------------------------------------------------------
    /** RawDictionaryEditorModule constructor. */

    public RawDictionaryEditorModule()
    {
        moduleHost_      = null;
        dictionaryMap_   = null;
        moduleID_        = -1;
        initialised_     = false;
        propertySupport_ = new PropertyChangeSupport(this);
    }

    //--------------------------------------------------------------------------

    /** Returns the name of this module. */
    public String toString()
    {
        return moduleName_;
    }

    //--------------------------------------------------------------------------

    public void addPropertyChangeListener(PropertyChangeListener listener)
    {
        propertySupport_.addPropertyChangeListener(listener);
    }

    //--------------------------------------------------------------------------

    public void removePropertyChangeListener(PropertyChangeListener listener)
    {
        propertySupport_.removePropertyChangeListener(listener);
    }

    //--------------------------------------------------------------------------
    //---- IModule Interface
    //--------------------------------------------------------------------------

    public int getModuleID()
    {
        return moduleID_;
    }

    //--------------------------------------------------------------------------

    public String getModuleName()
    {
        return moduleName_;
    }

    //--------------------------------------------------------------------------

    public String getModuleDescription()
    {
        return moduleDescription_;
    }

    //--------------------------------------------------------------------------

    public javax.swing.Action[] getActions()
    {
        return moduleActions_;
    }

    //--------------------------------------------------------------------------

    public boolean initialiseModule(IModuleHost moduleHost, int moduleID)
    {
        boolean bRet = false;

        try
        {
            // Prevent multiple initialisation.
            if (!initialised_)
            {
                // Store module Id and the reference to the module host.
                moduleHost_ = moduleHost;
                moduleID_   = moduleID;

                // Get a reference to the current case.
                ICaseServer currentCase = moduleHost_.getCaseServer();

                // Get list of dictionaries from the case's application class
                // object.
                String[] dictList =
                    currentCase.application().dictionaries();

                // Create map object to hold descriptor objects.
                dictionaryMap_ = new Hashtable();

                // Add Dictionaries node to tree.
                int dictRootNode =
                    moduleHost_.addNode
                    (
                        0,
                        "RawDictionaries",
                        FoamXTypeEx.getIcon(FoamXType.Type_Dictionary)
                    );

                // Create a sub-node for all dictionaries.
                for (int i = 0; i <dictList.length; i++)
                {
                    // Strip out the path for the display string.

                    String dictPath = dictList[i];
                    String dictName = dictPath;
                    int nIndex = dictPath.lastIndexOf("/");
                    if (nIndex != -1)
                    {
                        dictName = dictPath.substring(nIndex + 1);
                    }

                    // Get the DictionaryEntry reference for this dictionary.
                    IDictionaryEntryHolder holder = new IDictionaryEntryHolder();
                    currentCase.getDictionary(dictPath, holder);
                    IDictionaryEntry dictRoot = holder.value;

                    // construct the full path (eg 'system/controlDict')
                    dictPath =
                        dictRoot.typeDescriptor().dictionaryPath()
                        + "/" + dictPath;

                    // Create a descriptor object for this dictionary and
                    // add to map.
                    DictionaryDescriptor dictDesc =
                        new DictionaryDescriptor(dictName, dictPath);
                    dictionaryMap_.put(dictName, dictDesc);

                    // Add node for this dictionary.
                    int nodeID =
                        moduleHost_.addNode
                        (
                            dictRootNode,
                            dictName,
                            FoamXTypeEx.getIcon(FoamXType.Type_Dictionary)
                        );

                    // Attach the Open, Close and Save commands to this node.
                    moduleHost_.addNodeAction
                    (
                         nodeID,
                         moduleActions_[OPENRAWDICTIONARY_INDEX],
                         dictName
                    );
                    moduleHost_.addNodeAction
                    (
                        nodeID,
                        moduleActions_[CLOSERAWDICTIONARY_INDEX],
                        dictName
                    );
                    moduleHost_.addNodeAction
                    (
                        nodeID,
                        moduleActions_[SAVERAWDICTIONARY_INDEX],
                        dictName
                    );
                }

                // Indicate success.
                initialised_ = true;
                bRet          = true;
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return bRet;
    }

    //--------------------------------------------------------------------------

    public boolean shutdownModule()
    {
        return true;
    }

    //--------------------------------------------------------------------------
    //---- Visibility Interface
    //--------------------------------------------------------------------------

    public boolean avoidingGui()
    {
        return false;
    }

    //--------------------------------------------------------------------------

    public void dontUseGui()
    {
    }

    //--------------------------------------------------------------------------

    public boolean needsGui()
    {
        return true;
    }

    //--------------------------------------------------------------------------

    public void okToUseGui()
    {
    }

    //--------------------------------------------------------------------------

    protected void frameClosing(int nWindowID)
    {
        try
        {
            // Find the dictionary for this window.
            Enumeration iter = dictionaryMap_.elements();
            while (iter.hasMoreElements())
            {
                DictionaryDescriptor dictDesc =
                    (DictionaryDescriptor)iter.nextElement();
                if (dictDesc.getWindowID() == nWindowID)
                {
                    // See if the dictionary has been modified.
                    DictionaryTextEditor editor = dictDesc.getTextEditor();
                    if (editor.isModified())
                    {
                    }

                    // Reset dictionary descriptor.
                    dictDesc.setWindowID(-1);
                    dictDesc.setTextEditor(null);
                    break;
                }
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

    class OpenRawDictionaryAction
    extends AbstractAction
    {
        OpenRawDictionaryAction()
        {
            super(OPENRAWDICTIONARY_ACTION);
        }
        OpenRawDictionaryAction(String nm)
        {
            super(nm);
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                String dictName = e.getActionCommand();
                if (!dictionaryMap_.containsKey(dictName))
                {
                    throw new Exception
                    (
                        "Illegal dictionary name in OpenRawDictionary action."
                    );
                }

                // Get the dictionary descriptor object.
                DictionaryDescriptor dictDesc =
                    (DictionaryDescriptor)dictionaryMap_.get(dictName);

                // See if we already have a window for this dictionary.
                if (dictDesc.getWindowID() != -1)
                {
                    // Pop this window to the front.
                    moduleHost_.showWindow(dictDesc.getWindowID());
                }
                else
                {
                    // Need to create a new window for this dictionary.

                    // Get the dictionary contents from the current case.
                    org.omg.CORBA.StringHolder contentsHolder =
                        new org.omg.CORBA.StringHolder();
                    ICaseServer currentCase = moduleHost_.getCaseServer();
                    currentCase.readRawDictionary
                    (
                        dictDesc.getPath(),
                        contentsHolder
                    );

                    // Create a DictionaryEditor object.
                    DictionaryTextEditor dictEditor =
                        new DictionaryTextEditor(dictDesc.getName());
                    dictEditor.setEditorText(contentsHolder.value);

                    // Create a window to house the dictionary editor.
                    int nID = moduleHost_.createWindow(dictEditor);
                    if (nID != -1)
                    {
                        // Get the FoamXInternalFrame object.
                        FoamXInternalFrame frame =
                            moduleHost_.getWindowFrame(nID);
                        if (frame != null)
                        {
                            // Subscribe to the FoamXInternalFrame object's
                            // closing event.
                            frame.addInternalFrameListener
                            (
                                new javax.swing.event.InternalFrameAdapter()
                                {
                                    public void internalFrameClosing
                                    (
                                        javax.swing.event.InternalFrameEvent evt
                                    )
                                    {
                                        FoamXInternalFrame frame =
                                            (FoamXInternalFrame)evt.getSource();
                                        if (frame != null)
                                        {
                                            frameClosing(frame.getWindowID());
                                        }
                                    }
                                }
                            );

                            // Store reference to editor object and window ID.
                            dictDesc.setTextEditor(dictEditor);
                            dictDesc.setWindowID(nID);
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

    class SaveRawDictionaryAction
    extends AbstractAction
    {
        SaveRawDictionaryAction()
        {
            super(SAVERAWDICTIONARY_ACTION);
        }
        SaveRawDictionaryAction(String nm)
        {
            super(nm);
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                String dictName = e.getActionCommand();
                if (!dictionaryMap_.containsKey(dictName))
                {
                    throw new Exception
                    (
                        "Illegal dictionary name in SaveRawDictionary action."
                    );
                }

                App.printMessage("SaveRawDictionaryAction invoked!!");
                //moduleHost_.ShowModuleGUI(nModuleID_);
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    class CloseRawDictionaryAction
    extends AbstractAction
    {
        CloseRawDictionaryAction()
        {
            super(CLOSERAWDICTIONARY_ACTION);
        }
        CloseRawDictionaryAction(String nm)
        {
            super(nm);
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                String dictName = e.getActionCommand();
                if (!dictionaryMap_.containsKey(dictName))
                {
                    throw new Exception
                    (
                        "Illegal dictionary name in CloseRawDictionary action."
                    );
                }

                // Get the dictionary descriptor object.
                DictionaryDescriptor dictDesc =
                    (DictionaryDescriptor)dictionaryMap_.get(dictName);

                // See if we already have a window for this dictionary.
                if (dictDesc.getWindowID() != -1)
                {
                    // See if the dictionary has been modified.
                    DictionaryTextEditor editor = dictDesc.getTextEditor();
                    if (editor.isModified())
                    {
                    }

                    // Close the dictionary editor window. This will fire a
                    // Frame Closing event which
                    // we handle in FrameClosing where we tidy up the
                    // dictionary descriptor object.
                    moduleHost_.closeWindow(dictDesc.getWindowID());
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
}
