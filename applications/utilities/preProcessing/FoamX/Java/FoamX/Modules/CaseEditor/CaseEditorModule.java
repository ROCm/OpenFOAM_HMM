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
package FoamX.Modules.CaseEditor;

import java.util.*;
import java.net.*;
import java.io.*;
import java.beans.*;
import javax.swing.*;
import javax.swing.filechooser.*;
import java.awt.event.*;

import org.omg.CORBA.StringHolder;

import FoamXServer.*;
import FoamXServer.CaseServer.*;
import FoamXServer.CaseBrowser.*;

import FoamX.App;
import FoamX.Exceptions.FoamXException;
import FoamX.Util.FoamXTypeEx;
import FoamX.Util.FileEditor;
import FoamX.Util.FileMonitor;
import FoamX.Util.FileListener;
import FoamX.Util.FileEvent;
import FoamX.CaseManagement.IModule;
import FoamX.CaseManagement.IModuleHost;
import FoamX.CaseManagement.CaseChooserDlg;
import FoamX.CaseManagement.CaseBrowserPanel;
import FoamX.Reporting.ReportingWindow;

import FoamX.Editors.DictionaryEntryEditor.CompoundEntryListener;
import FoamX.Editors.DictionaryEntryEditor.CompoundEntryEvent;
import FoamX.Editors.DictionaryEntryEditor.DictionaryEditorPanel;

import FoamX.WindowManagement.FoamXInternalFrame;
import FoamX.Editors.ApplicationEditor.BoundaryDefinitionModelItem;
import FoamX.Editors.TypeEditor.TypeDescriptorCache;

public class CaseEditorModule
    extends Object
    implements
        IModule,                // Can be loaded as module
        PatchStatusListener,    // Gets warned of editing of patch types
        CompoundEntryListener,  // Gets warned of editing of dictionaries
        FileListener            // Gets warned if file on disk changes
{
    //--------------------------------------------------------------------------

    protected static final String MODULENAME = "CaseEditorModule";
    protected static final String MODULEDESCRIPTION =
        "Generic Case Editor Module.";

    protected PropertyChangeSupport propertySupport_;
    protected IModuleHost moduleHost_;
    protected int moduleID_;

    // Map of dictionary name to DictionaryWindowHandler objects.
    protected java.util.Hashtable dictionaryMap_;
    // Map of entry name to DictionaryEntryWindowHandler objects.
    protected java.util.Hashtable dictionaryEntryMap_;
    // Map of field name to FieldWindowHandler objects.
    protected java.util.Hashtable fieldWindowMap_;
    // Map of patch name to PatchWindowHandler objects.
    protected java.util.Hashtable patchWindowMap_;
    // Map of field and patch name to DictionaryWindowHandler objects.
    protected java.util.Hashtable patchfieldWindowMap_;

    // Listen for changes to dictionaries and fields
    protected FileMonitor monitor_;

    protected int meshNode_;
    protected int patchesNode_;

    //--------------------------------------------------------------------------
    // Module actions.
    javax.swing.Action[] moduleActions_;

    OpenDictionaryAction openDictionaryAction_;
    CloseDictionaryAction closeDictionaryAction_;
    SaveDictionaryAction saveDictionaryAction_;
    ReadMeshAction readMeshAction_;
    ImportMeshAction importMeshAction_;
    EditDictionaryEntryAction editDictionaryEntryAction_;
    EditPatchAction editPatchAction_;
    EditFieldAction editFieldAction_;
    EditPatchFieldParametersAction editPatchFieldParametersAction_;

    //--------------------------------------------------------------------------
    /** DictionaryEditorModule constructor. */
    public CaseEditorModule()
    {
        try
        {
            propertySupport_     = new PropertyChangeSupport(this);
            moduleHost_          = null;
            moduleID_            = -1;
            dictionaryMap_       = null;
            dictionaryEntryMap_  = null;
            fieldWindowMap_      = null;
            patchfieldWindowMap_ = null;
            monitor_             = null;
            meshNode_            = -1;
            patchesNode_         = -1;

            // Initialise the actions.
            openDictionaryAction_           = new OpenDictionaryAction();
            closeDictionaryAction_          = new CloseDictionaryAction();
            saveDictionaryAction_           = new SaveDictionaryAction();
            readMeshAction_                 = new ReadMeshAction();
            importMeshAction_               = new ImportMeshAction();
            editDictionaryEntryAction_      = new EditDictionaryEntryAction();
            editPatchAction_                = new EditPatchAction();
            editFieldAction_                = new EditFieldAction();
            editPatchFieldParametersAction_ =
                new EditPatchFieldParametersAction();

            moduleActions_ = new Action[8];
            moduleActions_[0] = openDictionaryAction_;
            moduleActions_[1] = closeDictionaryAction_;
            moduleActions_[2] = saveDictionaryAction_;
            moduleActions_[3] = readMeshAction_;
            moduleActions_[4] = importMeshAction_;
            moduleActions_[5] = editDictionaryEntryAction_;
            moduleActions_[6] = editPatchAction_;
            moduleActions_[7] = editFieldAction_;
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Returns the name of this module. */
    public String toString()
    {
        return MODULENAME;
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
        return MODULENAME;
    }

    //--------------------------------------------------------------------------

    public String getModuleDescription()
    {
        return MODULEDESCRIPTION;
    }

    //--------------------------------------------------------------------------

    public javax.swing.Action[] getActions()
    {
        return moduleActions_;
    }

    //--------------------------------------------------------------------------

    // Called by CasePanel to initialise.
    public boolean initialiseModule(IModuleHost moduleHost, int moduleID)
        throws FoamXIOError, FoamXError, Exception
    {
        boolean bRet = false;

        try
        {
            if (moduleID_ != -1)
            {
                throw new FoamXException("Multiple initialisation error.");
            }

            Icon fieldIcon   = App.getResources().getIcon("FieldImage");
            Icon meshIcon    = App.getResources().getIcon("MeshImage");
            Icon patchesIcon = App.getResources().getIcon("PatchesImage");

            // Store module Id and the reference to the module host.
            moduleHost_ = moduleHost;
            moduleID_   = moduleID;

            // Create map objects to hold window handler objects.
            dictionaryMap_       = new Hashtable();
            dictionaryEntryMap_  = new Hashtable();
            fieldWindowMap_      = new Hashtable();
            patchWindowMap_      = new Hashtable();
            patchfieldWindowMap_ = new Hashtable();

            // Get a temporary reference to the case server object.
            ICaseServer caseServer = moduleHost_.getCaseServer();

            //-----------------------------------------------------------------
            // Get list of dictionaries from the case's application class
            // object.
            String[] dictList = caseServer.application().dictionaries();

            // Add Dictionaries node to tree.
            int dictRootNode = moduleHost_.addNode
            (
                0,
                "Dictionaries",
                FoamXTypeEx.getIcon(FoamXType.Type_Dictionary)
            );

            // Create sub-nodes for all top-level dictionaries.
            for (int i = 0; i <dictList.length; i++)
            {
                String dictPath = dictList[i];

                // Create window handler for dictionary
                initialiseDictionaryWindowHandler(dictPath, dictRootNode);
            }

            //-----------------------------------------------------------------
            // Add Field nodes to tree.
            int fieldsRoot = moduleHost_.addNode(0, "Fields", fieldIcon);

            // Attach the Read command to this node.
            moduleHost_.addNodeAction
            (
                fieldsRoot,
                readMeshAction_,
                ""
            );

            String[] fieldNameList = caseServer.application().fields();
            IGeometricFieldHolder holder = new IGeometricFieldHolder();

            for (int nField = 0; nField <fieldNameList.length; nField++)
            {
                String fieldName = fieldNameList[nField];

                // Get field value object and add construct a
                // FieldWindowHandler object.
                caseServer.getFieldValues(fieldName, holder);

                // Create a window data object for this field and add to map.
                FieldWindowHandler windowData =
                    new FieldWindowHandler(fieldName, holder.value);
                fieldWindowMap_.put(fieldName, windowData);

                // Add node for this patch.
                int nodeID =
                    moduleHost_.addNode(fieldsRoot, fieldName, fieldIcon);

                // Attach the Edit command to this node.
                moduleHost_.addNodeAction
                (
                    nodeID,
                    editFieldAction_,
                    fieldName
                );
            }

            //-----------------------------------------------------------------
            // Add Mesh and Patch nodes to tree.
            meshNode_ = moduleHost_.addNode(0, "Mesh", meshIcon);
            patchesNode_ =
                moduleHost_.addNode
                (
                    meshNode_,
                    "Patches",
                    patchesIcon
                );

            // Attach the Mesh commands (read, import) to these nodes.
            moduleHost_.addNodeAction(meshNode_, readMeshAction_, "");
            moduleHost_.addNodeAction(patchesNode_, readMeshAction_, "");

            moduleHost_.addNodeAction(meshNode_, importMeshAction_, "");
            moduleHost_.addNodeAction(patchesNode_, importMeshAction_, "");

            // If there is a mesh defined, read it.
            if (caseServer.meshDefined())
            {
                readMeshAction_.actionPerformed(null);
            }

            //-----------------------------------------------------------------
            // Initialise the patchField window handlers.
            initialisePatchFieldWindowHandlers();

            // Indicate success.
            bRet = true;
        }
        catch (FoamXIOError ioErr)
        {
            App.printMessage
            (
                App.DEBUGLEVEL_VERBOSE,
                "CaseEditor:initialiseModule:FoamXIOError"
            );
            throw ioErr;
        }
        catch (FoamXError fxErr)
        {
            App.printMessage
            (
                App.DEBUGLEVEL_VERBOSE,
                "CaseEditor:initialiseModule:FoamXError"
            );
            throw fxErr;
        }
        catch (Exception ex)
        {
            App.printMessage
            (
                App.DEBUGLEVEL_VERBOSE,
                "CaseEditor:initialiseModule:Exception"
            );
            throw ex;
        }

        return bRet;
    }

    //--------------------------------------------------------------------------

    public boolean shutdownModule()
    {
        try
        {
            // Close all open windows managed by this module.

            // Close dictionary windows.
            Enumeration iter = dictionaryMap_.elements();
            while (iter.hasMoreElements())
            {
                DictionaryWindowHandler windowHandler =
                    (DictionaryWindowHandler)iter.nextElement();

                if (windowHandler.getWindowID() != -1)
                {
                    moduleHost_.closeWindow(windowHandler.getWindowID());
                    windowHandler.setWindowID(-1);
                }
            }

            // Close field windows.
            iter = fieldWindowMap_.elements();
            while (iter.hasMoreElements())
            {
                FieldWindowHandler windowHandler =
                    (FieldWindowHandler)iter.nextElement();

                if (windowHandler.getWindowID() != -1)
                {
                    moduleHost_.closeWindow(windowHandler.getWindowID());
                    windowHandler.setWindowID(-1);
                }
            }

            // Close patchField windows.
            iter = patchfieldWindowMap_.elements();
            while (iter.hasMoreElements())
            {
                DictionaryWindowHandler windowHandler =
                    (DictionaryWindowHandler)iter.nextElement();
                if (windowHandler.getWindowID() != -1)
                {
                    moduleHost_.closeWindow(windowHandler.getWindowID());
                    windowHandler.setWindowID(-1);
                }
            }

            // Close patch windows.
            iter = patchWindowMap_.elements();
            while (iter.hasMoreElements())
            {
                PatchWindowHandler windowHandler =
                    (PatchWindowHandler)iter.nextElement();

                if (windowHandler.getWindowID() != -1)
                {
                    moduleHost_.closeWindow(windowHandler.getWindowID());
                    windowHandler.setWindowID(-1);
                }
            }

            // Remove any file listeners
            unInstallFileListeners();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return true;
    }

    //--------------------------------------------------------------------------

    public boolean activate()
    {
        // Remove any file listeners
        installFileListeners();

        // Always return success
        return true;
    }
        
    //--------------------------------------------------------------------------

    public boolean deactivate()
    {
        // Remove any file listeners
        unInstallFileListeners();

        // Always return success
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

    public void okToUseGui()
    {
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
    //---- PropertyChangeListener Methods
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
    //---- PatchStatusListener Interface
    //--------------------------------------------------------------------------

    public void patchPhysicalTypeChanged(PatchStatusEvent evt)
    {
        try
        {
            // Get all open field windows to refresh themselves.
            {
                Enumeration iter = fieldWindowMap_.elements();
                while (iter.hasMoreElements())
                {
                    FieldWindowHandler windowHandler =
                        (FieldWindowHandler)iter.nextElement();
                    if (windowHandler.getWindowID() != -1)
                    {
                        //windowHandler.getWindow().refresh();
                        moduleHost_.closeWindow(windowHandler.getWindowID());
                        windowHandler.setWindowID(-1);
                    }
                }
            }

            // Close all open patchField windows.
            {
                Enumeration iter = patchfieldWindowMap_.elements();
                while (iter.hasMoreElements())
                {
                    DictionaryWindowHandler windowHandler =
                        (DictionaryWindowHandler)iter.nextElement();
                    if (windowHandler.getWindowID() != -1)
                    {
                        moduleHost_.closeWindow(windowHandler.getWindowID());
                        windowHandler.setWindowID(-1);
                    }
                }
            }

            // Re-initialise the patchField window handlers.
            initialisePatchFieldWindowHandlers();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //---- CompoundEntryListener Interface
    //--------------------------------------------------------------------------

    public boolean openSubDictionary(CompoundEntryEvent evt)
    {
        switch (evt.getType())
        {
            case CompoundEntryEvent.TYPE_DICTIONARY:
                {
                    // Open the sub-dictionary window.
                    ActionEvent evt2 =
                        new ActionEvent
                        (
                            this,
                            ActionEvent.ACTION_PERFORMED,
                            evt.getName()
                        );
                    openDictionaryAction_.actionPerformed(evt2);
                }
                break;
            case CompoundEntryEvent.TYPE_PATCH:
                {
                    // Open the patch parameters window.
                    ActionEvent evt2 =
                        new ActionEvent
                        (
                            this,
                            ActionEvent.ACTION_PERFORMED,
                            evt.getName()
                        );
                    editPatchFieldParametersAction_.actionPerformed(evt2);
                }
                break;
            default:
                App.handleException
                (
                    new Exception("Unknown event type in openSubDictionary.")
                );
        }

        return true;
    }

    //--------------------------------------------------------------------------
    //---- FileListener Interface
    //--------------------------------------------------------------------------

    public void fileOpened(FileEvent evt)
    {}

    public void fileClosed(FileEvent evt)
    {}

    public void fileChanged(FileEvent evt)
    {
        String fName = evt.file().getPath();

        String dictName = evt.file().getName();

        if (dictionaryMap_.containsKey(dictName))
        {
            int ret = JOptionPane.showConfirmDialog
            (
                App.getRootFrame(),
                "Dictionary " + dictName + " changed on disk\n"
              + "Do you want to reread it?",
                "FoamX...",
                JOptionPane.YES_NO_OPTION,
                JOptionPane.QUESTION_MESSAGE
            );
            if (ret == JOptionPane.OK_OPTION)
            {
                try
                {
                    // See if we already have a window for this dictionary
                    // (just to know whether we need to reopen one)
                    int oldWindowID = 
                        (
                            (DictionaryWindowHandler)dictionaryMap_.get
                            (
                                dictName
                            )
                        ).getWindowID();

                    // Close the dictionary window.
                    ActionEvent evt2 =
                        new ActionEvent
                        (
                            this,
                            ActionEvent.ACTION_PERFORMED,
                            dictName
                        );
                    closeDictionaryAction_.actionPerformed(evt2);

                    // Reload dictionary and recreate windowHandler
                    ICaseServer caseServer = moduleHost_.getCaseServer();

                    // Get the DictionaryEntry reference for this dictionary.
                    // (force rereading from disk)
                    IDictionaryEntryHolder holder =
                        new IDictionaryEntryHolder();
                    caseServer.getDictionary(dictName, true, holder);
                    IDictionaryEntry dictRoot = holder.value;

                    // Create a window handler object for this dictionary and
                    // add to map.
                    DictionaryWindowHandler windowHandler
                         = new DictionaryWindowHandler(dictName, dictRoot);

                    dictionaryMap_.put(dictName, windowHandler);

                    // Reopen window if it was open before
                    if (oldWindowID != -1)
                    {
                        openDictionaryAction_.actionPerformed(evt2);
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
        else
        {
            // Is not dictionary so is field or patch. Reload mesh to be on the
            // safe side. (no need to read mesh for change of field but alas)

            int ret = JOptionPane.showConfirmDialog
            (
                App.getRootFrame(),
                "File " + fName + " changed on disk\n"
              + "Do you want to reread Mesh&Fields?",
                "FoamX...",
                JOptionPane.YES_NO_OPTION,
                JOptionPane.QUESTION_MESSAGE
            );
            if (ret == JOptionPane.OK_OPTION)
            {
                readMeshAction_.actionPerformed(null);
            }
        }
    }

    //--------------------------------------------------------------------------
    //---- Window Handler Classes
    //--------------------------------------------------------------------------
    // Inner class for storing dictionary window information.
    private class DictionaryWindowHandler
    extends javax.swing.event.InternalFrameAdapter
    {
        private String dictionaryName_;
        private IDictionaryEntry dictRootEntry_;
        private int windowID_;
        private DictionaryEditorPanel windowObj_;

        public DictionaryWindowHandler
        (
            String dictionaryName,
            IDictionaryEntry rootEntry
        )
        {
            dictionaryName_ = dictionaryName;
            dictRootEntry_  = rootEntry;
            windowID_       = -1;
            windowObj_      = null;
        }

        public String getName()
        {
            return dictionaryName_;
        }

        public IDictionaryEntry getRootEntry()
        {
            return dictRootEntry_;
        }

        public int getWindowID()
        {
            return windowID_;
        }

        public void setWindowID(int windowID)
        {
            windowID_ = windowID;
        }

        public DictionaryEditorPanel getWindow()
        {
            return windowObj_;
        }

        public void setWindow(DictionaryEditorPanel windowObj)
        {
            windowObj_ = windowObj;
        }

        /**
         * Event handler for frame closing event. Resets the window ID and
         * window object reference.
         */
        public void internalFrameClosing
        (
            javax.swing.event.InternalFrameEvent evt
        )
        {
            FoamXInternalFrame frame = (FoamXInternalFrame)evt.getSource();
            if (frame != null)
            {
                windowID_  = -1;
                windowObj_ = null;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Inner class for storing dictionary entry window information.
    private class DictionaryEntryWindowHandler
    {
        private String dictionaryName_;
        private String itemName_;
        private IDictionaryEntry entry_;

        public DictionaryEntryWindowHandler
        (
            String dictionaryName,
            String itemName,
            IDictionaryEntry entry
        )
        {
            dictionaryName_ = dictionaryName;
            itemName_       = itemName;
            entry_          = entry;
        }

        public String getDictionaryName()
        {
            return dictionaryName_;
        }

        public String getItemName()
        {
            return itemName_;
        }

        public IDictionaryEntry getEntry()
        {
            return entry_;
        }
    }

    //--------------------------------------------------------------------------
    // Inner class for storing field window information.
    private class FieldWindowHandler
    extends javax.swing.event.InternalFrameAdapter
    {
        private String fieldName_;
        private IGeometricField fieldValues_;
        private int windowID_;
        private FieldPropertiesWindow windowObj_;

        public FieldWindowHandler
        (
            String fieldName,
            IGeometricField fieldValues
        )
        {
            fieldName_   = fieldName;
            fieldValues_ = fieldValues;
            windowID_    = -1;
            windowObj_   = null;
        }

        public String getName()
        {
            return fieldName_;
        }

        public IGeometricField getFieldValues()
        {
            return fieldValues_;
        }

        public int getWindowID()
        {
            return windowID_;
        }

        public void setWindowID(int windowID)
        {
            windowID_ = windowID;
        }

        public FieldPropertiesWindow getWindow()
        {
            return windowObj_;
        }

        public void setWindow(FieldPropertiesWindow windowObj)
        {
            windowObj_ = windowObj;
        }

        /**
         * Event handler for frame closing event. Resets the window ID and
         * window object reference.
         */
        public void internalFrameClosing
        (
            javax.swing.event.InternalFrameEvent evt
        )
        {
            FoamXInternalFrame frame = (FoamXInternalFrame)evt.getSource();
            if (frame != null)
            {
                windowID_  = -1;
                windowObj_ = null;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Inner class for storing patch window information.
    private class PatchWindowHandler
    extends javax.swing.event.InternalFrameAdapter
    {
        private String patchName_;
        private String patchPhysicalType_;
        private int windowID_;
        private PatchEditorWindow windowObj_;

        public PatchWindowHandler(String patchName, String patchPhysicalType)
        {
            patchName_    = patchName;
            patchPhysicalType_ = patchPhysicalType;
            windowID_     = -1;
            windowObj_    = null;
        }

        public String getName()
        {
            return patchName_;
        }

        public String getPatchPhysicalType()
        {
            return patchPhysicalType_;
        }

        public void setPatchPhysicalType(String patchPhysicalType)
        {
            patchPhysicalType_ = patchPhysicalType;
        }

        public int getWindowID()
        {
            return windowID_;
        }

        public void setWindowID(int windowID)
        {
            windowID_ = windowID;
        }

        public PatchEditorWindow getWindow()
        {
            return windowObj_;
        }

        public void setWindow(PatchEditorWindow windowObj)
        {
            windowObj_ = windowObj;
        }

        /**
         * Event handler for frame closing event. Resets the window ID and
         * window object reference.
         */
        public void internalFrameClosing
        (
            javax.swing.event.InternalFrameEvent evt
        )
        {
            FoamXInternalFrame frame = (FoamXInternalFrame)evt.getSource();
            if (frame != null)
            {
                windowID_  = -1;
                windowObj_ = null;
            }
        }
    }

    //--------------------------------------------------------------------------
    //---- Utility Methods
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // Install fileMonitor on all files this case uses
    private void installFileListeners()
    {
        try
        {
            ICaseServer caseServer = moduleHost_.getCaseServer();
            ICaseBrowser caseBrowser = moduleHost_.getCaseBrowser();

            //-----------------------------------------------------------------
            // Get list of dictionaries from the case's application class
            // object.
            String[] dictList = caseServer.application().dictionaries();
            String[] fieldNameList = caseServer.application().fields();

            // All dictionaries, field and mesh files
            File[] files = new File[dictList.length + fieldNameList.length + 1];

            int fileI = 0;

            String casePath =
                caseServer.caseRoot() + '/' + caseServer.caseName() + '/';

            for (int i = 0; i < dictList.length; i++)
            {
                // Get the DictionaryEntry reference for this dictionary.
                IDictionaryEntryHolder holder = new IDictionaryEntryHolder();
                caseServer.getDictionary(dictList[i], false, holder);
                IDictionaryEntry dictRoot = holder.value;

                String fName =
                    casePath
                  + dictRoot.typeDescriptor().dictionaryPath()
                  + '/'
                  + dictRoot.typeDescriptor().name();

                files[fileI++] = new File(fName);
            }

            //-----------------------------------------------------------------
            // Get field file names
            String time = caseServer.getTime();

            for (int i = 0; i < fieldNameList.length; i++)
            {
                String fName = casePath + time + '/' + fieldNameList[i];

                files[fileI++] = new File(fName);
            }

            //-----------------------------------------------------------------
            // Get mesh file names (only boundary file is used)
            files[fileI++] = new File(casePath + "constant/polyMesh/boundary");

            // Remove old file monitor
            unInstallFileListeners();

            // Install new file change listener
            int sleep = Integer.parseInt
            (
                App.getOptions().getProperty("FoamX.FileSleep", "1000")
            );

            monitor_ = new FileMonitor(caseBrowser, files, sleep);

            monitor_.addFileListener(this);

            Thread monitorThread = new Thread(monitor_);

            monitorThread.start();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    // Uninstall fileMonitor on all files this case uses
    private void unInstallFileListeners()
    {
        if (monitor_ != null)
        {
            monitor_.removeFileListener(this);
            monitor_.stop();
            monitor_ = null;
        }
    }

    //--------------------------------------------------------------------------

    private void printCaseMessage(String message)
    {
        // Get reporting window for this case.
        ICaseServer caseServer = moduleHost_.getCaseServer();
        String key =
            caseServer.caseRoot() + '/' + caseServer.caseName();
        ReportingWindow reportWin =
            App.getReportingManager().getReportingWindow(key);
        if (reportWin != null)
        {
            // Print message into main reporting window.
            reportWin.printMessage(message);
        }
    }

    //--------------------------------------------------------------------------

    private void AddChildren
    (
        int parentNodeID,
        String dictKey,
        IDictionaryEntry compoundEntry
    )
    {
        try
        {
            IDictionaryEntry[] subElements = compoundEntry.subElements();

            // Loop over all sub entries.
            for (int i=0; i <subElements.length; i++)
            {
                // Get TypeDescriptor for this entry.
                IDictionaryEntry entry          = subElements[i];
                ITypeDescriptor  typeDescriptor = entry.typeDescriptor();

                // Make sure that we have a name to display.
                String entryName = typeDescriptor.displayName();
                if (entryName == null || entryName.length() == 0)
                {
                    entryName = typeDescriptor.name();
                }

                // Add node for this entry.
                FoamXType type = typeDescriptor.type();
                int nodeID = moduleHost_.addNode
                (
                    parentNodeID,
                    entryName,
                    TypeDescriptorCache.getIcon(typeDescriptor)
                );

                // If this entry is another dictionary, add its children to
                // the tree.
                if (type == FoamXType.Type_Dictionary)
                {
                    // Get the DictionaryEntry reference for this dictionary.
                    IDictionaryEntry dictEntry =
                        IDictionaryEntryHelper.narrow(entry);

                    // Create a descriptor object for this dictionary and add
                    // to list.
                    String subDictName = dictEntry.typeDescriptor().name();

                    if (subDictName.length()> 0)
                    {
                        // Sub dictionary key names are constructed from the
                        // dictionary names delimited with "/".
                        String key = dictKey + "/" + subDictName;

                        // Attach the Open and Close commands to this node.
                        moduleHost_.addNodeAction
                        (
                            nodeID,
                            openDictionaryAction_,
                            key
                        );
                        moduleHost_.addNodeAction
                        (
                            nodeID,
                            closeDictionaryAction_,
                            key
                        );

                        // Create a window handler for this dictionary and add
                        // to map.
                        DictionaryWindowHandler windowHandler =
                            new DictionaryWindowHandler
                            (
                                subDictName,
                                dictEntry
                            );
                        dictionaryMap_.put(key, windowHandler);

                        // Recursively add this dictionary's children to the
                        // tree.
                        AddChildren(nodeID, key, dictEntry);
                    }
                }
                else
                {
                    // Dictionary entry key names are constructed from the
                    // dictionary name and the item name delimited with "::".
                    String key = dictKey + "::" + entryName;

                    // Attach the Edit command to this node.
                    moduleHost_.addNodeAction
                    (
                        nodeID,
                        editDictionaryEntryAction_,
                        key
                    );

                    IDictionaryEntry ncEntry =
                        IDictionaryEntryHelper.narrow(entry);

                    // Create a window handler for this dictionary entry and
                    // add to map.
                    DictionaryEntryWindowHandler windowHandler =
                        new DictionaryEntryWindowHandler
                        (
                            dictKey,
                            entryName,
                            ncEntry
                        );
                    dictionaryEntryMap_.put(key, windowHandler);
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    // Create window handlers for a dictionary and its subdictionaries.
    private void initialiseDictionaryWindowHandler
    (
        String dictPath,
        int dictRootNode
    )
    {
        try
        {
            // Strip out the directory to get the name.
            String dictName = dictPath;
            int nIndex = dictPath.lastIndexOf("/");
            if (nIndex != -1)
            {
                dictName = dictPath.substring(nIndex + 1);
            }

            ICaseServer caseServer = moduleHost_.getCaseServer();

            // Get the DictionaryEntry reference for this dictionary.
            IDictionaryEntryHolder holder = new IDictionaryEntryHolder();
            caseServer.getDictionary(dictPath, false, holder);
            IDictionaryEntry dictRoot = holder.value;

            // Create a window handler object for this dictionary and add
            // to map.
            DictionaryWindowHandler windowHandler
                 = new DictionaryWindowHandler(dictName, dictRoot);
            dictionaryMap_.put(dictName, windowHandler);

            // Add node for this dictionary.
            int nodeID = moduleHost_.addNode
            (
                dictRootNode,
                dictName,
                TypeDescriptorCache.getIcon(dictRoot.typeDescriptor())
            );

            // Attach the Open, Close, Save commands to this
            // node.
            moduleHost_.addNodeAction
            (
                nodeID,
                openDictionaryAction_,
                dictName
            );
            moduleHost_.addNodeAction
            (
                nodeID,
                closeDictionaryAction_,
                dictName
            );
            moduleHost_.addNodeAction
            (
                nodeID,
                saveDictionaryAction_,
                dictName
            );

            // Recursively add this dictionary's children to the tree.
            AddChildren(nodeID, dictName, dictRoot);
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

    //--------------------------------------------------------------------------

    // Create window handlers for all patches of all fields.
    private void initialisePatchFieldWindowHandlers()
    {
        try
        {
            // Clear all patch field parameter data.
            patchfieldWindowMap_.clear();

            // Get a temporary reference to the case server object.
            ICaseServer caseServer = moduleHost_.getCaseServer();

            // Loop over all fields and check whether parameters are required
            // for this patch.
            String[] fieldNameList = caseServer.application().fields();
            for (int nField = 0; nField <fieldNameList.length; nField++)
            {
                String fieldName = fieldNameList[nField];

                // Get field values object for the specified field.
                IGeometricFieldHolder holder = new IGeometricFieldHolder();
                caseServer.getFieldValues(fieldName, holder);

                IGeometricField fieldValObject = holder.value;

                IDictionaryEntryHolder dictEntryHolder =
                    new IDictionaryEntryHolder();

                // Loop over all currently defined patches.
                String[] patchNames = caseServer.patchNames();
                for (int nPatch = 0; nPatch <patchNames.length; nPatch++)
                {
                    String patchName = patchNames[nPatch];

                    // Get dictionary entry object for the extra parameters
                    // if required.
                    fieldValObject.getPatchFieldParameters
                    (
                        patchName,
                        dictEntryHolder
                    );
                    if (dictEntryHolder.value != null)
                    {
                        // Use fieldName/patchName as key.
                        String key = fieldName + "/" + patchName;

                        // Create a window handler for this dictionary and add
                        // to map.
                        DictionaryWindowHandler windowHandler =
                            new DictionaryWindowHandler
                            (
                                key,
                                dictEntryHolder.value
                            );
                        patchfieldWindowMap_.put(key, windowHandler);

                        // Recursively add this dictionary's children to the
                        // tree.
                        //AddChildren(nodeID, key, dictEntry);
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

    private void initialisePatchWindowHandlers()
    {
        try
        {
            // Close patch windows.
            Enumeration iter = patchWindowMap_.elements();
            while (iter.hasMoreElements())
            {
                PatchWindowHandler windowHandler =
                    (PatchWindowHandler)iter.nextElement();

                if (windowHandler.getWindowID() != -1)
                {
                    moduleHost_.closeWindow(windowHandler.getWindowID());
                    windowHandler.setWindowID(-1);
                }
            }
            patchWindowMap_.clear();


            // Add patch nodes (patchnames from caseServer)

            Icon patchIcon   = App.getResources().getIcon("PatchImage");

            ICaseServer caseServer = moduleHost_.getCaseServer();

            if (caseServer.meshDefined())
            {
                String[] patchNameList = caseServer.patchNames();
                for
                (
                    int nPatch = 0;
                    nPatch < patchNameList.length;
                    nPatch++
                )
                {
                    String patchName = patchNameList[nPatch];

                    // Get patch boundary type.
                    StringHolder patchPhysicalTypeHolder = new StringHolder();
                    caseServer.getPatchPhysicalType
                    (
                        patchName,
                        patchPhysicalTypeHolder
                    );

                    // Create a window data object for this patch and add
                    // to map.
                    PatchWindowHandler windowData =
                        new PatchWindowHandler
                        (
                            patchName,
                            patchPhysicalTypeHolder.value
                        );
                    patchWindowMap_.put(patchName, windowData);

                    // Add node for this patch.
                    int nodeID =
                        moduleHost_.addNode
                        (
                            patchesNode_,
                            patchName,
                            patchIcon
                        );

                    // Attach the Edit command to this node.
                    moduleHost_.addNodeAction
                    (
                        nodeID,
                        editPatchAction_,
                        patchName
                    );
                }
            }
            else
            {
                JOptionPane.showMessageDialog
                (
                    App.getRootFrame(),
                    "No mesh defined.",
                    "FoamX...",
                    JOptionPane.WARNING_MESSAGE
                );
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

    class OpenDictionaryAction
    extends AbstractAction
    {
        OpenDictionaryAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("OpenDictionaryImage")
            );
            putValue(Action.NAME, "Open Dictionary");
            putValue(Action.SHORT_DESCRIPTION, "Open Dictionary");
            putValue(Action.LONG_DESCRIPTION, "Open Dictionary");
        }

        public synchronized void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get dictionary name stored in action command string.
                String dictName = e.getActionCommand();

                if (!dictionaryMap_.containsKey(dictName))
                {
                    throw new Exception
                    (
                        "Illegal dictionary name in OpenDictionary action."
                    );
                }

                // Get the dictionary window handler object.
                DictionaryWindowHandler windowHandler =
                    (DictionaryWindowHandler)dictionaryMap_.get(dictName);

                // See if we already have a window for this dictionary.
                if (windowHandler.getWindowID() != -1)
                {
                    // Pop this window to the front.
                    moduleHost_.showWindow(windowHandler.getWindowID());
                }
                else
                {
                    // Need to create a new DictionaryEditorWindow for this
                    // dictionary.

                    // Get the dictionary's root IDictionaryEntry interface.
                    IDictionaryEntry rootEntry = windowHandler.getRootEntry();

                    // Create a DictionaryEditorPanel object.
                    // Need to pass the key to this dictionary so that the
                    // user can open sub-item windows from this window.
                    DictionaryEditorPanel dictPanel =
                        new DictionaryEditorPanel(rootEntry, dictName);

                    // Create a window to house the dictionary editor.
                    int windowID = moduleHost_.createWindow(dictPanel);
                    if (windowID != -1)
                    {
                        // Print a message
                        printCaseMessage
                        (
                            "Opening dictionary " + dictName +
                            " at time " + moduleHost_.getCaseServer().getTime()
                        );

                        // Get the FoamXInternalFrame object.
                        FoamXInternalFrame frame =
                            moduleHost_.getWindowFrame(windowID);
                        if (frame != null)
                        {
                            // Store window ID and window reference.
                            windowHandler.setWindowID(windowID);
                            windowHandler.setWindow(dictPanel);

                            // Subscribe to the FoamXInternalFrame object's
                            // closing event.
                            // Direct events to the window handler object.
                            frame.addInternalFrameListener(windowHandler);

                            // Subscribe to the dictionary editor panel's
                            // openSubDictionary event.
                            // Direct events to the module.
                            dictPanel.
                                getTableModel().
                                    addCompoundEntryListener
                                    (
                                        CaseEditorModule.this
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
    }

    //--------------------------------------------------------------------------

    class CloseDictionaryAction
    extends AbstractAction
    {
        CloseDictionaryAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("CloseDictionaryImage")
            );
            putValue(Action.NAME, "Close Dictionary");
            putValue(Action.SHORT_DESCRIPTION, "Close Dictionary");
            putValue(Action.LONG_DESCRIPTION, "Close Dictionary");
        }

        public synchronized void actionPerformed(ActionEvent e)
        {
            try
            {
                String dictName = e.getActionCommand();
                if (!dictionaryMap_.containsKey(dictName))
                {
                    throw new Exception
                    (
                        "Illegal dictionary name in CloseDictionary action."
                    );
                }

                // Get the dictionary window handler object.
                DictionaryWindowHandler windowHandler =
                    (DictionaryWindowHandler)dictionaryMap_.get(dictName);

                // See if we already have a window for this dictionary.
                if (windowHandler.getWindowID() != -1)
                {
                    // Close the dictionary editor window. This will fire a
                    // Frame Closing event which
                    // we handle in the window handler where we tidy up.
                    moduleHost_.closeWindow(windowHandler.getWindowID());
                    windowHandler.setWindowID(-1);
                }
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    class SaveDictionaryAction
    extends AbstractAction
    {
        SaveDictionaryAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("SaveDictionaryImage")
            );
            putValue(Action.NAME, "Save Dictionary");
            putValue(Action.SHORT_DESCRIPTION, "Save Dictionary");
            putValue(Action.LONG_DESCRIPTION, "Save Dictionary");
        }

        public synchronized void actionPerformed(ActionEvent e)
        {
            try
            {
                String dictName = e.getActionCommand();
                if (!dictionaryMap_.containsKey(dictName))
                {
                    throw new Exception
                    (
                        "Illegal dictionary name in Save Dictionary action."
                    );
                }

                // Get the dictionary window handler object.
                DictionaryWindowHandler windowHandler =
                    (DictionaryWindowHandler)dictionaryMap_.get(dictName);

                // Get the dictionary's root IDictionaryEntry interface.
                IDictionaryEntry rootEntry = windowHandler.getRootEntry();

                // Remove any listeners (should be for this dictionary only
                // but alas)
                unInstallFileListeners();

                // Save the root dictionary.
                rootEntry.save();

                // Reinstall listeners
                installFileListeners();
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

    // Start file editor. File passed in event. Not (yet) used.
    class EditRawDictionaryAction
    extends AbstractAction
    {
        EditRawDictionaryAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("EditRawDictionaryImage")
            );
            putValue(Action.NAME, "Edit Dictionary (raw)");
            putValue(Action.SHORT_DESCRIPTION, "Edit Dictionary Raw");
            putValue(Action.LONG_DESCRIPTION, "Edit Dictionary Raw");
        }

        public synchronized void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get dictionary name stored in action command string.
                String fileName = e.getActionCommand();

                new FileEditor
                (
                    false,      // non modal
                    false,      // edit mode
                    null,       // no caseBrowser available
                    fileName,
                    0,
                    0,
                    0
                );
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    class ReadMeshAction
    extends AbstractAction
    {
        ReadMeshAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("ReadMeshImage")
            );
            putValue(Action.NAME, "Read Mesh&Fields");
            putValue(Action.SHORT_DESCRIPTION, "Read Mesh&Fields");
            putValue(Action.LONG_DESCRIPTION, "Read Mesh&Fields");
        }

        public synchronized void actionPerformed(ActionEvent e)
        {
            try
            {
                Icon patchesIcon = App.getResources().getIcon("PatchesImage");

                // Tidy up any previous mesh data. Since we don't have
                // access to pathWindowIds themselves delete and reconstruct
                // whole patch node
                moduleHost_.removeNode(patchesNode_);
                patchesNode_ =
                    moduleHost_.addNode
                    (
                        meshNode_,
                        "Patches",
                        patchesIcon
                    );

                // Get a temporary reference to the case server object.
                ICaseServer caseServer = moduleHost_.getCaseServer();

                if (caseServer.meshDefined())
                {
                    // Get the case server to reread the mesh and fields
                    caseServer.readMesh();

                    // Might have changed the number of patches
                    initialisePatchWindowHandlers();

                    patchPhysicalTypeChanged
                    (
                        // Dummy patch status event
                        new PatchStatusEvent
                        (
                            this,
                            "ReadMeshAction",
                            "ReadMeshAction"
                        )
                    );
                }
                else
                {
                    JOptionPane.showMessageDialog
                    (
                        App.getRootFrame(),
                        "No mesh defined.",
                        "FoamX...",
                        JOptionPane.WARNING_MESSAGE
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

    class ImportMeshAction
    extends AbstractAction
    {
        ImportMeshAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("ReadMeshImage")
            );
            putValue(Action.NAME, "Import Mesh");
            putValue(Action.SHORT_DESCRIPTION, "Import Mesh");
            putValue(Action.LONG_DESCRIPTION, "Import Mesh");
        }

        public synchronized void actionPerformed(ActionEvent e)
        {
            try
            {
                CaseChooserDlg caseChooser =
                    new CaseChooserDlg
                    (
                        App.getRootFrame(),
                        CaseBrowserPanel.SELECT_CASE_MODE
                    );

                caseChooser.show();

                if 
                (
                    (caseChooser.getCaseRoot() == null)
                 || (caseChooser.getCaseName() == null)
                 || (caseChooser.getCaseBrowser() == null)
                )
                {
                    return;
                }
 
                String rootAndCase =
                     caseChooser.getCaseRoot()
                     + "/" + caseChooser.getCaseName();

                if
                (
                    JOptionPane.showConfirmDialog
                    (
                        App.getRootFrame(),
                        "Are you sure you want to import mesh from\n"
                        + rootAndCase + "?"
                        + "\n(and overwrite current mesh)",
                        "Confirm Import Operation",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.QUESTION_MESSAGE
                    ) != JOptionPane.OK_OPTION
                )
                {
                    return;
                }

                // Get a temporary reference to the case server object.
                ICaseServer caseServer = moduleHost_.getCaseServer();

                System.out.println
                (
                    "Import mesh from" + rootAndCase
                     + " on machine " + caseChooser.getHostName()
                );

                caseServer.importMesh
                (
                    caseChooser.getHostName(),
                    caseChooser.getCaseRoot(),
                    caseChooser.getCaseName()
                );

                // Force mesh reread
                readMeshAction_.actionPerformed(e);
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

//    class ImportMeshDescriptionAction
//    extends AbstractAction
//    {
//        ImportMeshDescriptionAction()
//        {
//            putValue
//            (
//                Action.SMALL_ICON,
//                App.getResources().getIcon("ReadMeshImage")
//            );
//            putValue(Action.NAME, "Import meshDescription");
//            putValue(Action.SHORT_DESCRIPTION, "Import meshDescription");
//            putValue(Action.LONG_DESCRIPTION, "Import meshDescription");
//        }
//
//        public synchronized void actionPerformed(ActionEvent e)
//        {
//            try
//            {
//                CaseChooserDlg caseChooser =
//                    new CaseChooserDlg
//                  (
//                        App.getRootFrame(),
//                        CaseBrowserPanel.SELECT_CASE_MODE
//                    );
//                caseChooser.show();
//                if 
//                (
//                    (caseChooser.getCaseRoot() == null)
//                 || (caseChooser.getCaseName() == null)
//                 || (caseChooser.getCaseBrowser() == null)
//                )
//                {
//                    return;
//                }
//
//                String rootAndCase =
//                    caseChooser.getCaseRoot()
//                    + "/" + caseChooser.getCaseName();
//
//                if
//                (
//                    JOptionPane.showConfirmDialog
//                    (
//                        App.getRootFrame(),
//                        "Are you sure you want to import meshDescription"
//                        + " file from\n"
//                        + rootAndCase + "?"
//                        + "\n(and overwrite any current mesh)",
//                        "Confirm Import Operation",
//                        JOptionPane.YES_NO_OPTION,
//                        JOptionPane.QUESTION_MESSAGE
//                    ) != JOptionPane.OK_OPTION
//                )
//                {
//                    return;
//                }
//
//                StringHolder holder = new StringHolder();
//                caseChooser.getCaseBrowser().getHostname(holder);
//                String hostName = holder.value;
//
//                System.out.println("Copy " + hostName + ":" + rootAndCase);
//
////                //Hack: copying through Corba link for now.
////                StringHolder holder = new StringHolder();
////                caseChooser.getCaseBrowser().readFile
////                (
////                    rootAndCase + "/" + "constant/polyMesh/meshDescription",
////                    holder
////                );
////                String contents = holder.value;
////
////                // Get a temporary reference to the case server object.
////                ICaseServer caseServer = moduleHost_.getCaseServer();
////
////                caseServer.writeFile
////                (
////                    "constant/polyMesh/meshDescription",
////                    contents
////                );
//
//                // Get a temporary reference to the case server object.
//                ICaseServer caseServer = moduleHost_.getCaseServer();
//                caseServer.importMeshDescription
//                (
//                    caseChooser.getHostName(),
//                    caseChooser.getCaseRoot(),
//                    caseChooser.getCaseName()
//                );
//
//                // Force mesh reread
//                readMeshAction_.actionPerformed(e);
//            }
//            catch (Exception ex)
//            {
//                App.handleAllExceptions(ex);
//            }
//        }
//    }
//
    //--------------------------------------------------------------------------

    class EditDictionaryEntryAction
    extends AbstractAction
    {
        EditDictionaryEntryAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("EditDictionaryEntryImage")
            );
            putValue(Action.NAME, "Edit Dictionary Entry");
            putValue(Action.SHORT_DESCRIPTION, "Edit Dictionary Entry");
            putValue(Action.LONG_DESCRIPTION, "Edit Dictionary Entry");
        }

        public synchronized void actionPerformed(ActionEvent e)
        {
            try
            {
                String entryName = e.getActionCommand();
                if (!dictionaryEntryMap_.containsKey(entryName))
                {
                    throw new Exception
                    (
                        "Illegal entry name in EditDictionaryEntry action."
                    );
                }

                // Get the dictionary entry window handler object.
                DictionaryEntryWindowHandler entryHandler =
                    (DictionaryEntryWindowHandler)dictionaryEntryMap_.get
                        (entryName);

                // Get the parent dictionary window handler object.
                String dictName = entryHandler.getDictionaryName();
                if (!dictionaryMap_.containsKey(dictName))
                {
                    throw new Exception
                    (
                        "Illegal dictionary name in EditDictionaryEntry action."
                    );
                }
                DictionaryWindowHandler dictHandler =
                    (DictionaryWindowHandler)dictionaryMap_.get(dictName);

                // See if we already have a window for this dictionary.
                if (dictHandler.getWindowID() != -1)
                {
                    // Pop this window to the front.
                    moduleHost_.showWindow(dictHandler.getWindowID());

                    // Select the appropriate dictionary entry.
                    dictHandler.getWindow().selectEntry
                    (
                        entryHandler.getItemName()
                    );
                }
                else
                {
                    // Need to create a new DictionaryEditorWindow for this
                    // dictionary.

                    // Invoke the Open Dictionary action.
                    ActionEvent evt2 =
                        new ActionEvent
                        (
                            this,
                            ActionEvent.ACTION_PERFORMED,
                            dictName
                        );
                    openDictionaryAction_.actionPerformed(evt2);

                    // Should have a window open now, so select the
                    // appropriate entry
                    if (dictHandler.getWindowID() != -1)
                    {
                        // Select the appropriate dictionary entry.
                        dictHandler.getWindow().selectEntry
                        (
                            entryHandler.getItemName()
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

    class EditPatchAction
    extends AbstractAction
    {
        EditPatchAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("EditPatchImage")
            );
            putValue(Action.NAME, "Define Boundary Type");
            putValue(Action.SHORT_DESCRIPTION, "Define Boundary Type");
            putValue(Action.LONG_DESCRIPTION, "Define Boundary Type");
        }

        public synchronized void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get patch name stored in action command string.
                String patchName = e.getActionCommand();
                if (!patchWindowMap_.containsKey(patchName))
                {
                    throw new Exception
                    (
                        "Illegal patch name in EditPatchAction action."
                    );
                }

                // Get the patch window handler object.
                PatchWindowHandler windowHandler =
                    (PatchWindowHandler)patchWindowMap_.get(patchName);

                // See if we already have a window for this patch.
                if (windowHandler.getWindowID() != -1)
                {
                    // Pop this window to the front.
                    moduleHost_.showWindow(windowHandler.getWindowID());
                }
                else
                {
                    // Need to create a new PatchEditorWindow for this patch.

                    // Get a reference to the current case.
                    ICaseServer caseServer = moduleHost_.getCaseServer();

                    // Create a patch editor window object.
                    PatchEditorWindow patchEditorWindow =
                        new PatchEditorWindow(caseServer, patchName);

                    // Create a window to house the field value editor.
                    int windowID =
                        moduleHost_.createWindow(patchEditorWindow);
                    if (windowID != -1)
                    {
                        // Get the FoamXInternalFrame object.
                        FoamXInternalFrame frame =
                            moduleHost_.getWindowFrame(windowID);
                        if (frame != null)
                        {
                            // Store window ID and window reference.
                            windowHandler.setWindowID(windowID);
                            windowHandler.setWindow(patchEditorWindow);

                            // Subscribe to the FoamXInternalFrame object's
                            // closing event.
                            // Direct events to the window handler object.
                            frame.addInternalFrameListener(windowHandler);

                            // Listen out for patch status changes.
                            // Direct events to the module.
                            patchEditorWindow.addPatchStatusListener
                            (
                                CaseEditorModule.this
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
    }

    //--------------------------------------------------------------------------

    class EditFieldAction
    extends AbstractAction
    {
        EditFieldAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("EditFieldImage")
            );
            putValue(Action.NAME, "Edit Field");
            putValue(Action.SHORT_DESCRIPTION, "Edit Field");
            putValue(Action.LONG_DESCRIPTION, "Edit Field");
        }

        public synchronized void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get field name stored in action command string.
                String fieldName = e.getActionCommand();
                if (!fieldWindowMap_.containsKey(fieldName))
                {
                    throw new Exception
                    (
                        "Illegal field name in EditFieldAction action."
                    );
                }

                // Get the field window handler object.
                FieldWindowHandler windowHandler =
                    (FieldWindowHandler)fieldWindowMap_.get(fieldName);

                // See if we already have a window for this field.
                if (windowHandler.getWindowID() != -1)
                {
                    // Pop this window to the front.
                    moduleHost_.showWindow(windowHandler.getWindowID());
                }
                else
                {
                    // Get a reference to the current case.
                    ICaseServer caseServer = moduleHost_.getCaseServer();

                    // Create a new field window.
                    FieldPropertiesWindow fieldWindow =
                        new FieldPropertiesWindow(caseServer, fieldName);

                    // Create a window to house the field value editor.
                    int windowID = moduleHost_.createWindow(fieldWindow);
                    if (windowID != -1)
                    {
                        // Print a message
                        printCaseMessage
                        (
                            "Opening field " + caseServer.getTime() + '/'
                             + fieldName
                        );

                        // Get the FoamXInternalFrame object.
                        FoamXInternalFrame frame =
                            moduleHost_.getWindowFrame(windowID);
                        if (frame != null)
                        {
                            // Store window ID and window reference.
                            windowHandler.setWindowID(windowID);
                            windowHandler.setWindow(fieldWindow);

                            // Subscribe to the FoamXInternalFrame object's
                            // closing event.
                            // Direct events to the window handler object.
                            frame.addInternalFrameListener(windowHandler);

                            // Subscribe to the field window's
                            // openSubDictionary event.
                            // Direct events to the module.
                            fieldWindow.getTableModel()
                                .addCompoundEntryListener
                                (
                                    CaseEditorModule.this
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
    }

    //--------------------------------------------------------------------------

    class EditPatchFieldParametersAction
    extends AbstractAction
    {
        EditPatchFieldParametersAction()
        {
            putValue
            (
                Action.SMALL_ICON,
                App.getResources().getIcon("EditFieldImage")
            );
            putValue(Action.NAME, "Edit Patch Field Parameters");
            putValue(Action.SHORT_DESCRIPTION, "Edit Patch Field Parameters");
            putValue(Action.LONG_DESCRIPTION, "Edit Patch Field Parameters");
        }

        public synchronized void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get field and patch names stored in action command string.
                String patchFieldName = e.getActionCommand();

                if (!patchfieldWindowMap_.containsKey(patchFieldName))
                {
                    throw new Exception
                    (
                        "Illegal patch/field name in " +
                        "EditPatchFieldParametersAction action."
                    );
                }

                // Get the dictionary window handler object.
                DictionaryWindowHandler windowHandler =
                    (DictionaryWindowHandler)patchfieldWindowMap_.get
                    (
                        patchFieldName
                    );

                // See if we already have a window for this dictionary.
                if (windowHandler.getWindowID() != -1)
                {
                    // Pop this window to the front.
                    moduleHost_.showWindow(windowHandler.getWindowID());
                }
                else
                {
                    // Need to create a new DictionaryEditorWindow for this
                    // dictionary.

                    // Get the dictionary's root IDictionaryEntry interface.
                    IDictionaryEntry rootEntry = windowHandler.getRootEntry();

                    // Create a DictionaryEditorPanel object.
                    // Need to pass the key to this dictionary so that the
                    // user can open sub-item windows from this window.
                    DictionaryEditorPanel dictPanel =
                        new DictionaryEditorPanel(rootEntry, patchFieldName);

                    // Create a window to house the dictionary editor.
                    int windowID = moduleHost_.createWindow(dictPanel);
                    if (windowID != -1)
                    {
                        // Get the FoamXInternalFrame object.
                        FoamXInternalFrame frame =
                            moduleHost_.getWindowFrame(windowID);
                        if (frame != null)
                        {
                            // Store window ID and window reference.
                            windowHandler.setWindowID(windowID);
                            windowHandler.setWindow(dictPanel);

                            // Subscribe to the FoamXInternalFrame object's
                            // closing event.
                            // Direct events to the window handler object.
                            frame.addInternalFrameListener(windowHandler);

                            // Subscribe to the dictionary editor panel's
                            // openSubDictionary event.
                            // Direct events to the module.
                            dictPanel.getTableModel().addCompoundEntryListener
                            (
                                CaseEditorModule.this
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
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}







