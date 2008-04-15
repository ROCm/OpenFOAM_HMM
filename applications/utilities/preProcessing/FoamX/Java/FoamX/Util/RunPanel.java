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
package FoamX.Util;

import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import org.omg.CORBA.StringHolder;

import FoamXServer.*;
import FoamXServer.CaseBrowser.*;

import FoamX.App;
import FoamX.Exceptions.FoamXCancel;
import FoamX.Exceptions.FoamXException;
import FoamX.Util.BusyCursor;
import FoamX.Util.FileEditor;
import FoamX.Tools.InvokeUtilityAction;
import FoamX.CaseManagement.HostChooserDlg;
import FoamX.CaseManagement.CaseManager;
import FoamX.Util.FileMonitor;
import FoamX.Util.FileListener;
import FoamX.Util.FileEvent;

import FoamX.Editors.CompoundEditor;
import FoamX.Editors.DictionaryEntryEditor.EntryCache.DictionaryCache;
import FoamX.Editors.DictionaryEntryEditor.EntryCache.DictionaryEntryCache;

public class RunPanel
    extends javax.swing.JDialog
    implements FileListener
{
    //--------------------------------------------------------------------------

    private static final int DEFAULT_HEIGHT = 400;
    private static final int DEFAULT_WIDTH = 400;

    //--------------------------------------------------------------------------

    protected java.awt.Frame parent_;
    protected ICaseBrowser caseBrowser_;
    protected String caseRoot_;
    protected String caseName_;
    protected String appName_;
    protected String dictName_;
    protected DictionaryCache controlDict_;         // Utility top level dict
    protected boolean hasArgs_;                     // has 'arguments' entry
    protected boolean nonDefaultArgs_;              // needs other than 
                                                    // just case/root

    protected String[] args_;


    //--------------------------------------------------------------------------
    /** Creates new form RunPanel */
    public RunPanel
    (
        java.awt.Frame parent,
        ICaseBrowser caseBrowser,
        String caseRoot,
        String caseName,
        String appName,                     // application to invoke
        boolean isApp,                      // application or utility
        String dictName,                    // dictionary to edit (if any)
        DictionaryCache controlDict         // entry of dictionary (if any)
    )
    {
        super(parent, "Run Editor", false); // Non modal.

        parent_ = parent;

        try
        {
            caseBrowser_ = caseBrowser;
            caseRoot_    = caseRoot;
            caseName_    = caseName;
            appName_     = appName;
            dictName_    = dictName;
            controlDict_ = controlDict;

            initComponents();


            // Title
            setTitle(appName_);

            //
            // Initialise various input fields
            //

            // Enable background running for Applications
            if (isApp)
            {
                backgroundBox_.setSelected(true);
            }

            // Start without args
            argListText_.setText("");
            argListText_.setEnabled(false);

            // Run on local host
            StringHolder holder = new StringHolder();
            caseBrowser_.getHostName(holder);
            machineText_.setText(holder.value);

            args_ = new String[1];
            args_[0] = "";

            if (controlDict_ == null)
            {
                dictNameText_.setText("");
                dictNameText_.setEnabled(false);

                editDictButton_.setEnabled(false);

                // Set argument list to root and case
                setArgsField(true, null);
                editArgsButton_.setEnabled(false);
            }
            else
            {
                dictNameText_.setText(controlDict_.getEntryName());
                dictNameText_.setEnabled(false);

                editDictButton_.setEnabled(true);

                // Set argument list from entries in controlDict
                setDefaultArgs(controlDict_);
                //// Disalllow editing
                // editArgsButton_.setEnabled(nonDefaultArgs_);
            }

            if (isApp)
            {
                String rootAndCase = caseRoot_ + "/" + caseName_;
                String logName = rootAndCase + '/' + appName_ + ".log";
                logText_.setText(logName);
            }
            else
            {
                logText_.setText("");
            }


            show();
        }
        catch(FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch(Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    // Fill in root and case in arguments in controlDict and check if more
    // than root and case are needed.
    private void setDefaultArgs(DictionaryCache dict)
    {
        hasArgs_ = false;
        nonDefaultArgs_ = false;

        // Get the arguments subelement of the utility dictionary
        DictionaryEntryCache argumentsEntry = dict.getSubEntry("arguments");

        if (argumentsEntry != null)
        {
            hasArgs_ = true;
            nonDefaultArgs_ = true;

            // Loop over all argument entries and replace rootAndCase with
            // current values.
            DictionaryEntryCache[] args = argumentsEntry.getSubElements();

            for (int i = 0; i < args.length; i++)
            {
                DictionaryEntryCache arg = args[i];

                int typeValue = arg.getTypeDescriptor().getType().value();

                if (typeValue == FoamXType._Type_RootDir)
                {
                    arg.updateValue(caseRoot_);
                }
                else if (typeValue == FoamXType._Type_RootAndCase)
                {
                    arg.updateValue(caseRoot_ + "/" + caseName_);

                    if (args.length == 1)
                    {
                        // Only one argument which is case and root.
                        nonDefaultArgs_ = false;
                    }
                }
                else if (typeValue == FoamXType._Type_CaseName)
                {
                    arg.updateValue(caseName_);
                }
            }

            // Update displayed argument string
            setArgsField(false, argumentsEntry);
        }
    }

    //--------------------------------------------------------------------------
    // Set arguments textField from control dictionary arguments or
    // straight from root,case. 
    private void setArgsField
    (
        boolean isApp,
        DictionaryEntryCache argumentsEntry
    )
    {
        if (argumentsEntry != null)
        {
            //- Get arguments as vector first
            Vector argsVector = new Vector();
            argumentsEntry.toStringRaw(argsVector);

            // Set both args_(for calling) and argString(for displaying)

            args_ = new String[argsVector.size()];
            String argString = "";

            int argi = 0;

            for (int i = 0; i < argsVector.size(); i++)
            {
                String arg = (String)argsVector.elementAt(i);

                args_[argi++] = arg;
                argString += arg;

                if (i <= argsVector.size() - 1)
                {
                    argString += ' ';
                }
            }

            argListText_.setText(argString);
        }
        else if (isApp)
        {
            args_ = new String[2];
            args_[0] = "-case";
            args_[1] = caseRoot_ + '/' + caseName_;
            argListText_.setText(args_[0] + " " + args_[1]);
        }
    }

    //--------------------------------------------------------------------------
    // Has dictionary been modified (and needs to be saved)?
    // If only 'arguments' entry changed don't bother saving
    private boolean dictModified()
    {
        try
        {
            if (controlDict_ == null)
            {
                return false;
            }

            // Underlying dictionary has not been modified
            if (!controlDict_.getDictEntry().modified())
            {
                return false;
            }

            // Has been modified. Check if only arguments have been modified.

            int nSubEntries = controlDict_.getDictEntry().nSubElements();
            if (nSubEntries == 0)
            {
                return false;
            }
            if (nSubEntries > 1)
            {
                return true;
            }

            // Only one entry in controlDict. Is arguments?
            if (hasArgs_)
            {
                return false;
            }
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }

        return true;
    }

    //--------------------------------------------------------------------------
    //---- FileListener Interface
    //--------------------------------------------------------------------------

    public void fileOpened(FileEvent evt)
    {
    }

    public void fileClosed(FileEvent evt)
    {
    }

    public void fileChanged(FileEvent evt)
    {
        String fName = evt.file().getPath();

        int index = fName.lastIndexOf('.');

        if (index >= 0)
        {
            JOptionPane.showMessageDialog
            (
                this,
                "Application " + appName_
              + " (pid " + fName.substring(index+1) + ") finished",
                "FoamX...",
                JOptionPane.INFORMATION_MESSAGE
            );
        }
        FileMonitor fm = (FileMonitor)evt.getSource();

        fm.removeFileListener(this);
    }

    //--------------------------------------------------------------------------
    //---- GUI 
    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents()//GEN-BEGIN:initComponents
    {
        mainPanel_ = new javax.swing.JPanel();
        editDictButton_ = new javax.swing.JButton();
        dictNameText_ = new javax.swing.JTextField();
        editArgsButton_ = new javax.swing.JButton();
        argListText_ = new javax.swing.JTextField();
        startRunButton_ = new javax.swing.JButton();
        machineLabel_ = new javax.swing.JLabel();
        machineText_ = new javax.swing.JTextField();
        backgroundBox_ = new javax.swing.JCheckBox();
        logLabel_ = new javax.swing.JLabel();
        viewLogButton_ = new javax.swing.JButton();
        logText_ = new javax.swing.JTextField();
        buttonPanel_ = new javax.swing.JPanel();
        closeButton_ = new javax.swing.JButton();
        
        getContentPane().setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Run Editor");
        setFont(new java.awt.Font("Dialog", 0, 10));
        setName("Run Editor");
        addWindowListener(new java.awt.event.WindowAdapter()
        {
            public void windowClosing(java.awt.event.WindowEvent evt)
            {
                closeDialog(evt);
            }
        });
        
        mainPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        mainPanel_.setBorder(new javax.swing.border.EmptyBorder(new java.awt.Insets(10, 10, 10, 10)));
        mainPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        mainPanel_.setPreferredSize(new java.awt.Dimension(400, 200));
        editDictButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        editDictButton_.setText("Edit Dictionary");
        editDictButton_.setMaximumSize(new java.awt.Dimension(117, 24));
        editDictButton_.setMinimumSize(new java.awt.Dimension(117, 24));
        editDictButton_.setPreferredSize(new java.awt.Dimension(117, 21));
        editDictButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                editDictButton_ActionPerformed(evt);
            }
        });
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        mainPanel_.add(editDictButton_, gridBagConstraints2);
        
        dictNameText_.setFont(new java.awt.Font("Dialog", 0, 10));
        dictNameText_.setText("dictionaryName");
        dictNameText_.setMinimumSize(new java.awt.Dimension(250, 17));
        dictNameText_.setPreferredSize(new java.awt.Dimension(300, 17));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.gridwidth = 5;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        mainPanel_.add(dictNameText_, gridBagConstraints2);
        
        editArgsButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        editArgsButton_.setText("Edit Arguments");
        editArgsButton_.setMaximumSize(new java.awt.Dimension(117, 24));
        editArgsButton_.setMinimumSize(new java.awt.Dimension(117, 24));
        editArgsButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                editArgsButton_ActionPerformed(evt);
            }
        });
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 1;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(0, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        mainPanel_.add(editArgsButton_, gridBagConstraints2);
        
        argListText_.setFont(new java.awt.Font("Dialog", 0, 10));
        argListText_.setText("root case");
        argListText_.setMinimumSize(new java.awt.Dimension(250, 17));
        argListText_.setPreferredSize(new java.awt.Dimension(300, 17));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 1;
        gridBagConstraints2.gridwidth = 5;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        mainPanel_.add(argListText_, gridBagConstraints2);
        
        startRunButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        startRunButton_.setText("Execute");
        startRunButton_.setMaximumSize(new java.awt.Dimension(117, 24));
        startRunButton_.setMinimumSize(new java.awt.Dimension(117, 24));
        startRunButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                startRunButton_ActionPerformed(evt);
            }
        });
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 2;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(0, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        mainPanel_.add(startRunButton_, gridBagConstraints2);
        
        machineLabel_.setFont(new java.awt.Font("Dialog", 0, 10));
        machineLabel_.setText("machine:");
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 2;
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.EAST;
        mainPanel_.add(machineLabel_, gridBagConstraints2);
        
        machineText_.setFont(new java.awt.Font("Dialog", 0, 10));
        machineText_.setText("penfold");
        machineText_.setToolTipText("Machine Name");
        machineText_.setMinimumSize(new java.awt.Dimension(60, 15));
        machineText_.setPreferredSize(new java.awt.Dimension(60, 15));
        machineText_.addMouseListener(new java.awt.event.MouseAdapter()
        {
            public void mouseClicked(java.awt.event.MouseEvent evt)
            {
                machineText_MouseClicked(evt);
            }
        });
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 2;
        gridBagConstraints2.gridy = 2;
        gridBagConstraints2.gridwidth = 3;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        mainPanel_.add(machineText_, gridBagConstraints2);
        
        backgroundBox_.setFont(new java.awt.Font("Dialog", 0, 10));
        backgroundBox_.setText("background");
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 5;
        gridBagConstraints2.gridy = 2;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.EAST;
        mainPanel_.add(backgroundBox_, gridBagConstraints2);
        
        logLabel_.setFont(new java.awt.Font("Dialog", 0, 10));
        logLabel_.setText("log:");
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 3;
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.EAST;
        mainPanel_.add(logLabel_, gridBagConstraints2);
        
        viewLogButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        viewLogButton_.setText("View Log");
        viewLogButton_.setMaximumSize(new java.awt.Dimension(117, 24));
        viewLogButton_.setMinimumSize(new java.awt.Dimension(117, 24));
        viewLogButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                viewLogButton_ActionPerformed(evt);
            }
        });
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 3;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(0, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        mainPanel_.add(viewLogButton_, gridBagConstraints2);
        
        logText_.setFont(new java.awt.Font("Dialog", 0, 10));
        logText_.setText("log.icoFoam");
        logText_.setMinimumSize(new java.awt.Dimension(250, 15));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 2;
        gridBagConstraints2.gridy = 3;
        gridBagConstraints2.gridwidth = 4;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        mainPanel_.add(logText_, gridBagConstraints2);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        getContentPane().add(mainPanel_, gridBagConstraints1);
        
        buttonPanel_.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.RIGHT));
        
        buttonPanel_.setBorder(new javax.swing.border.EtchedBorder());
        buttonPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        buttonPanel_.setPreferredSize(new java.awt.Dimension(400, 100));
        closeButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        closeButton_.setText("Close");
        closeButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnClose(evt);
            }
        });
        
        buttonPanel_.add(closeButton_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints1.weightx = 1.0;
        getContentPane().add(buttonPanel_, gridBagConstraints1);
        
        pack();
        java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        setSize(new java.awt.Dimension(480, 200));
        setLocation((screenSize.width-480)/2,(screenSize.height-200)/2);
    }//GEN-END:initComponents

    private void machineText_MouseClicked(java.awt.event.MouseEvent evt)//GEN-FIRST:event_machineText_MouseClicked
    {//GEN-HEADEREND:event_machineText_MouseClicked
        // Add your handling code here:

        // Pop up host chooser panel
        HostChooserDlg hostChooser =
            new HostChooserDlg(parent_);
        hostChooser.show();

        if (hostChooser.getHostName() != null)
        {
            machineText_.setText(hostChooser.getHostName());
        }
    }//GEN-LAST:event_machineText_MouseClicked

    private void viewLogButton_ActionPerformed(java.awt.event.ActionEvent evt)//GEN-FIRST:event_viewLogButton_ActionPerformed
    {//GEN-HEADEREND:event_viewLogButton_ActionPerformed
        // Add your handling code here:

        if (logText_.getText().length() > 0)
        {
            new FileEditor
            (
                true,           // modal
                false,          // Edit mode
                caseBrowser_,
                logText_.getText(),
                0,
                0,
                0
            );
        }

    }//GEN-LAST:event_viewLogButton_ActionPerformed

    private void startRunButton_ActionPerformed(java.awt.event.ActionEvent evt)//GEN-FIRST:event_startRunButton_ActionPerformed
    {//GEN-HEADEREND:event_startRunButton_ActionPerformed
        // Add your handling code here:

        try
        {
            // Try to save first
            CaseManager caseManager = App.getCaseManager();
            if (caseManager == null)
            {
                throw new FoamXException
                (
                    "No caseManager in top-level App."
                );
            }
            caseManager.saveCaseNice(caseRoot_, caseName_);

            // Actually invoke utility on the caseBrowser.
            int pid = caseBrowser_.invokeUtility
            (
                machineText_.getText(),
                appName_,
                args_,
                logText_.getText(),
                backgroundBox_.isSelected()
            );

            if (backgroundBox_.isSelected())
            {
                // When running in background print pid started with and
                // install filemonitor to check for finishing.
                if (pid != -1)
                {
                    // Do like JOptionPane.showMessageDialog but non-modal.
                    // (from sable-vm)
                    JOptionPane pane =
                        new JOptionPane
                        (
                            "Application " + appName_
                          + " started with pid " + pid,
                            JOptionPane.INFORMATION_MESSAGE
                        );
                    JDialog dialog = pane.createDialog(this, "FoamX...");
                    dialog.setModal(false);
                    dialog.pack();
                    dialog.show();

                    // Install new file change listener
                    StringHolder holder = new StringHolder();
                    caseBrowser_.getEnv("FOAM_JOB_DIR", holder);

                    File[] files = new File[1];
                    files[0] =
                        new File
                        (
                            holder.value
                          + "/finishedJobs/"
                          + machineText_.getText()
                          + "."
                          + pid
                        );

                    int sleep = Integer.parseInt
                    (
                        App.getOptions().getProperty("FoamX.FileSleep", "5000")
                    );

                    FileMonitor finishedMonitor =
                        new FileMonitor(caseBrowser_, files, sleep);

                    finishedMonitor.addFileListener(this);

                    Thread monitorThread = new Thread(finishedMonitor);

                    monitorThread.start();
                }
            }
            else
            {
                // Foreground running. Warn when finished.
                JOptionPane.showMessageDialog
                (
                    this,
                    "Application " + appName_ + " (pid " + pid
                  + ") finished successfully",
                    "FoamX...",
                    JOptionPane.INFORMATION_MESSAGE
                );
            }
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch (FoamXException ex)
        {
            App.handleAllExceptions(ex);
        }

    }//GEN-LAST:event_startRunButton_ActionPerformed

    private void editArgsButton_ActionPerformed(java.awt.event.ActionEvent evt)//GEN-FIRST:event_editArgsButton_ActionPerformed
    {//GEN-HEADEREND:event_editArgsButton_ActionPerformed
        // Add your handling code here:

        try
        {
            if (controlDict_ != null)
            {
                // Get the arguments subelement of the utility dictionary
                DictionaryEntryCache argumentsEntry =
                    controlDict_.getSubEntry("arguments");
                if (argumentsEntry != null)
                {
                    CompoundEditor compoundEditor =
                        new CompoundEditor
                        (
                            parent_,
                            argumentsEntry.getDictEntry()
                        );
                    compoundEditor.show();

                    // Update displayed value
                    setArgsField(false, argumentsEntry);

                }
                if (dictModified())
                {
                    controlDict_.getDictEntry().save();
                }
            }
        }
        catch (FoamXIOError ioErr)
        {
            App.handleAllExceptions(ioErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }

    }//GEN-LAST:event_editArgsButton_ActionPerformed

    private void editDictButton_ActionPerformed(java.awt.event.ActionEvent evt)//GEN-FIRST:event_editDictButton_ActionPerformed
    {//GEN-HEADEREND:event_editDictButton_ActionPerformed
        // Add your handling code here:

        try
        {
            CompoundEditor compoundEditor =
                new CompoundEditor  
                (
                    parent_,
                    controlDict_.getDictEntry()
                );

            compoundEditor.show();

            DictionaryEntryCache argumentsEntry =
                controlDict_.getSubEntry("arguments");

            // Update displayed value
            setArgsField(false, argumentsEntry);

            if (dictModified())
            {
                controlDict_.getDictEntry().save();
            }
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch (FoamXIOError ioErr)
        {
            App.handleAllExceptions(ioErr);
        }

    }//GEN-LAST:event_editDictButton_ActionPerformed

    //--------------------------------------------------------------------------

  private void shapeTypeComboActionPerformed (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_shapeTypeComboActionPerformed
  }//GEN-LAST:event_shapeTypeComboActionPerformed

    //--------------------------------------------------------------------------

  private void OnClose (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnClose
        setVisible(false);
        dispose();
  }//GEN-LAST:event_OnClose

    //--------------------------------------------------------------------------
    /** Closes the dialog */

  private void closeDialog(java.awt.event.WindowEvent evt)
    {//GEN-FIRST:event_closeDialog
        setVisible(false);
        dispose();
  }//GEN-LAST:event_closeDialog

    //--------------------------------------------------------------------------

  // Variables declaration - do not modify//GEN-BEGIN:variables
  private javax.swing.JPanel mainPanel_;
  private javax.swing.JButton editDictButton_;
  private javax.swing.JTextField dictNameText_;
  private javax.swing.JButton editArgsButton_;
  private javax.swing.JTextField argListText_;
  private javax.swing.JButton startRunButton_;
  private javax.swing.JLabel machineLabel_;
  private javax.swing.JTextField machineText_;
  private javax.swing.JCheckBox backgroundBox_;
  private javax.swing.JLabel logLabel_;
  private javax.swing.JButton viewLogButton_;
  private javax.swing.JTextField logText_;
  private javax.swing.JPanel buttonPanel_;
  private javax.swing.JButton closeButton_;
  // End of variables declaration//GEN-END:variables

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



