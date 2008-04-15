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
package FoamX.ProcessManagement;


import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Vector;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;

import org.omg.CORBA.StringHolder;

import FoamX.App;
import FoamX.Editors.StringListEditor;
import FoamX.Util.TableSorter;

import FoamXServer.*;
import FoamXServer.CaseBrowser.ICaseBrowser;

public class ProcessEditor extends javax.swing.JDialog
{
    //--------------------------------------------------------------------------

    protected java.awt.Frame parent_;
    protected ICaseBrowser caseBrowser_;
    protected RunningProcessModel runningModel_;
    protected FinishedProcessModel finishedModel_;
    protected String licenceDir_;
    
    protected TableSorter runningSorter_;
    protected TableSorter finishedSorter_;

    //--------------------------------------------------------------------------
    /** Creates new Process Editor */
    public ProcessEditor(java.awt.Frame parent, ICaseBrowser caseBrowser)
    {
        super(parent, "Process Editor", false);   // non Modal.

        try
        {
            parent_ = parent;

            // Initialise the process Manager
            caseBrowser_ = caseBrowser;

            // Initialise the table models and associated sorters
            runningModel_ = new RunningProcessModel(caseBrowser_);
            finishedModel_ = new FinishedProcessModel(caseBrowser_);
            runningSorter_ = new TableSorter(runningModel_);
            finishedSorter_ = new TableSorter(finishedModel_);

            initComponents();


            // Make sorters lister to tableheader clicks
            runningSorter_.addMouseListenerToHeaderInTable(runningTable_);
            finishedSorter_.addMouseListenerToHeaderInTable(finishedTable_);

            // Set licence dir
            StringHolder holder = new StringHolder();
            caseBrowser.getEnv("FOAM_JOB_DIR", holder);
            licenceDir_ = holder.value;

            // Start off with disabled job-wise buttons
            runInfoButton_.setEnabled(false);
            endButton_.setEnabled(false);
            endNowButton_.setEnabled(false);
            killButton_.setEnabled(false);
            suspButton_.setEnabled(false);
            contButton_.setEnabled(false);
            finInfoButton_.setEnabled(false);
            removeButton_.setEnabled(false);

            // Set status of representation checkboxes
            myRunsBox_.setSelected(true);
            myRunsBox_ActionPerformed
            (
                new java.awt.event.ActionEvent
                (
                    this,
                    java.awt.event.ActionEvent.ACTION_PERFORMED,
                    ""
                )
            );
            compactBox_.setSelected(false);

            // Initialise models
            caseBrowser_.refreshJobsLists();
            runningModel_.refresh
            (
                myRunsBox_.isSelected(),
                compactBox_.isSelected()
            );
            finishedModel_.refresh
            (
                myRunsBox_.isSelected(),
                compactBox_.isSelected()
            );

            // Title
            setLicLabel();

            // Set initial size of individual columns.
            initColumnSizes(runningTable_, runningModel_);
            initTableSize(runningTable_, runningModel_);
            initColumnSizes(finishedTable_, finishedModel_);
            initTableSize(finishedTable_, finishedModel_);

            tabbedPanel_.setSelectedComponent(finPanel_);
            
            repaint();

    //        tabbedPanel_.setSelectedComponent(runPanel_);
    //        repaint();

    //        runningTable_.getColumn
    //        (
    //            runningModel_.getColumnName(ProcessModel.STAT_INDEX)
    //        ).setCellRenderer(new StatusCellRenderer());
    //        finishedTable_.getColumn
    //        (
    //            finishedModel_.getColumnName(ProcessModel.STAT_INDEX)
    //        ).setCellRenderer(new StatusCellRenderer());


            // Listen out for selection events and enable the appropriate
            // buttons.

            runningTable_.getSelectionModel().addListSelectionListener
            (
                new ListSelectionListener()
                {
                    public void valueChanged(ListSelectionEvent evt)
                    {
                        enableRunButtons();
                    }
                }
            );

            finishedTable_.getSelectionModel().addListSelectionListener
            (
                new ListSelectionListener()
                {
                    public void valueChanged(ListSelectionEvent evt)
                    {
                        enableFinButtons();
                    }
                }
            );
        }
        catch(FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch(FoamXIOError ioErr)
        {
            App.handleAllExceptions(ioErr);
        }
    }

    //--------------------------------------------------------------------------
    /** Confirmation dialog helper */
    private boolean askConfirm(String mesg)
    {
        return
            JOptionPane.showConfirmDialog
            (
                this,
                mesg,
                "Confirmation Dialog",
                JOptionPane.YES_NO_OPTION,
                JOptionPane.QUESTION_MESSAGE
            ) ==
            JOptionPane.OK_OPTION;
    }
    
    //--------------------------------------------------------------------------
    /** Set licLabel textField to show number of licenses */
    private void setLicLabel()
    {
        licenceLabel_.setText
        (
            "jobControl : " + licenceDir_
          + " (" + runningModel_.getNJobs() + " jobs"
          + ", " + runningModel_.getNCountedProcs() + " floatlicensed jobs)"
        );
    }

    
    //--------------------------------------------------------------------------
    /** Size columns. Uses longValues() method from table */

    static private void initColumnSizes(JTable table, ProcessModel model)
    {
        TableColumn column = null;
        java.awt.Component comp = null;
        int headerWidth = 0;
        int cellWidth = 0;
        Object[] longValues = model.getLongValues();

        for (int i = 0; i < model.getColumnCount(); i++)
        {
            column = table.getColumnModel().getColumn(i);

            // comp = column.getDefaultRenderer();
            // headerWidth = comp.getPreferredSize().width;
            comp =
                table.getDefaultRenderer(model.getColumnClass(i)).
                    getTableCellRendererComponent
                    (
                         table, longValues[i],
                         false, false, 0, i
                    );
            cellWidth = comp.getPreferredSize().width;

            // column.setPreferredWidth(Math.max(headerWidth, cellWidth));
            column.setPreferredWidth(cellWidth);
        }
    } 

    //--------------------------------------------------------------------------
    /** Size table itself */

    static private void initTableSize(JTable table, ProcessModel model)
    {
        int height = table.getRowHeight() * model.getRowCount();

        table.setPreferredSize
        (
            new java.awt.Dimension(table.getWidth(), height)
        );
        table.revalidate();
    } 

    //--------------------------------------------------------------------------
    /** selectionhandler for runningTable */
    private void enableRunButtons()
    {
        boolean validSelection =
            !runningTable_.getSelectionModel().isSelectionEmpty();

        runInfoButton_.setEnabled(validSelection);
        endNowButton_.setEnabled(validSelection);
        endButton_.setEnabled(validSelection);
        killButton_.setEnabled(validSelection);
        suspButton_.setEnabled(validSelection);
        contButton_.setEnabled(validSelection);
    }
    
    //--------------------------------------------------------------------------
    /** selectionhandler for finishedTable */
    private void enableFinButtons()
    {
        boolean validSelection =
            !finishedTable_.getSelectionModel().isSelectionEmpty();
        finInfoButton_.setEnabled(validSelection);
        removeButton_.setEnabled(validSelection);
    }

//    //--------------------------------------------------------------------------
//    /** Clicked on entry. Detect double clicks. */
//    private void onJobClicked(MouseEvent e)
//    {
//        if (e.getClickCount() == 2)
//        {
//            JTable table = (JTable)e.getComponent();
//            Point origin = e.getPoint();
//            int row = table.rowAtPoint(origin);
//            if (row != -1)
//            {
//                TableSorter sortedModel = (TableSorter)table.getModel();
//                ProcessModel model = (ProcessModel)sortedModel.getModel();
//
//                displayJobInfo(model, row);
//            }
//        }
//    }


    //--------------------------------------------------------------------------
    /** Clicked on entry. Detect double clicks. */
    private void displayJobInfo(ProcessModel model, int row)
    {
        int size = model.getColumnCount();
        Vector labels = new Vector(size);
        Vector values = new Vector(size);
        for (int i = 0; i < size; i++)
        {
            labels.addElement(model.getColumnName(i));
            values.addElement(model.getValueAt(row, i));
        }

        StringListEditor dlg =
            new StringListEditor
            (
                parent_,
                "Job Info",
                labels, values
            );

        dlg.show();
    }


    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents()//GEN-BEGIN:initComponents
    {
        headerPanel_ = new javax.swing.JPanel();
        licenceLabel_ = new javax.swing.JLabel();
        tabbedPanel_ = new javax.swing.JTabbedPane();
        runPanel_ = new javax.swing.JPanel();
        runAllPanel_ = new java.awt.Panel();
        readButton_ = new javax.swing.JButton();
        statusButton_ = new javax.swing.JButton();
        purgeButton_ = new javax.swing.JButton();
        runPidPanel_ = new java.awt.Panel();
        runInfoButton_ = new javax.swing.JButton();
        endNowButton_ = new javax.swing.JButton();
        endButton_ = new javax.swing.JButton();
        killButton_ = new javax.swing.JButton();
        suspButton_ = new javax.swing.JButton();
        contButton_ = new javax.swing.JButton();
        runTableScrollPane_ = new javax.swing.JScrollPane();
        runningTable_ = new javax.swing.JTable();
        finPanel_ = new javax.swing.JPanel();
        finAllPanel_ = new java.awt.Panel();
        purgeOldButton_ = new javax.swing.JButton();
        finPidPanel_ = new java.awt.Panel();
        finInfoButton_ = new javax.swing.JButton();
        removeButton_ = new javax.swing.JButton();
        finTableScrollPane_ = new javax.swing.JScrollPane();
        finishedTable_ = new javax.swing.JTable();
        closeButtonPanel_ = new javax.swing.JPanel();
        myRunsBox_ = new javax.swing.JCheckBox();
        compactBox_ = new javax.swing.JCheckBox();
        closeButton_ = new javax.swing.JButton();
        
        getContentPane().setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        setTitle("Process Editor");
        setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
        setName("Process Editor");
        setFont(new java.awt.Font("Dialog", 0, 10));
        addWindowListener(new java.awt.event.WindowAdapter()
        {
            public void windowClosing(java.awt.event.WindowEvent evt)
            {
                closeDialog(evt);
            }
        });
        
        headerPanel_.setPreferredSize(new java.awt.Dimension(172, 23));
        headerPanel_.setFont(new java.awt.Font("Dialog", 1, 12));
        licenceLabel_.setText("/aa/bb/cc/foam/licences");
        licenceLabel_.setForeground(java.awt.Color.black);
        licenceLabel_.setFont(new java.awt.Font("Dialog", 0, 12));
        licenceLabel_.setOpaque(true);
        headerPanel_.add(licenceLabel_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.gridwidth = 2;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        getContentPane().add(headerPanel_, gridBagConstraints1);
        
        tabbedPanel_.setPreferredSize(new java.awt.Dimension(735, 286));
        runPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        runPanel_.setBorder(new javax.swing.border.EmptyBorder(new java.awt.Insets(10, 10, 10, 10)));
        runPanel_.setToolTipText("Running Jobs");
        runPanel_.setPreferredSize(new java.awt.Dimension(730, 161));
        runPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        runAllPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints3;
        
        runAllPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        runAllPanel_.setName("All Processes");
        runAllPanel_.setBackground(new java.awt.Color(204, 204, 204));
        runAllPanel_.setForeground(java.awt.Color.black);
        readButton_.setToolTipText("Reread licence directory");
        readButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        readButton_.setText("read");
        readButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnRead(evt);
            }
        });
        
        gridBagConstraints3 = new java.awt.GridBagConstraints();
        gridBagConstraints3.gridx = 0;
        gridBagConstraints3.gridy = 0;
        gridBagConstraints3.fill = java.awt.GridBagConstraints.HORIZONTAL;
        runAllPanel_.add(readButton_, gridBagConstraints3);
        
        statusButton_.setToolTipText("Contact machines to update status");
        statusButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        statusButton_.setText("status");
        statusButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnStatus(evt);
            }
        });
        
        gridBagConstraints3 = new java.awt.GridBagConstraints();
        gridBagConstraints3.gridx = 1;
        gridBagConstraints3.gridy = 0;
        gridBagConstraints3.fill = java.awt.GridBagConstraints.HORIZONTAL;
        runAllPanel_.add(statusButton_, gridBagConstraints3);
        
        purgeButton_.setToolTipText("Remove non-running jobs");
        purgeButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        purgeButton_.setText("purge");
        purgeButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnPurge(evt);
            }
        });
        
        gridBagConstraints3 = new java.awt.GridBagConstraints();
        gridBagConstraints3.gridx = 2;
        gridBagConstraints3.gridy = 0;
        gridBagConstraints3.fill = java.awt.GridBagConstraints.HORIZONTAL;
        runAllPanel_.add(purgeButton_, gridBagConstraints3);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        runPanel_.add(runAllPanel_, gridBagConstraints2);
        
        runPidPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints4;
        
        runPidPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        runPidPanel_.setName("Single Process");
        runPidPanel_.setBackground(new java.awt.Color(204, 204, 204));
        runPidPanel_.setForeground(java.awt.Color.black);
        runInfoButton_.setToolTipText("Display job info");
        runInfoButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        runInfoButton_.setText("Info");
        runInfoButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                onRunInfo(evt);
            }
        });
        
        gridBagConstraints4 = new java.awt.GridBagConstraints();
        gridBagConstraints4.gridx = 0;
        gridBagConstraints4.gridy = 0;
        gridBagConstraints4.fill = java.awt.GridBagConstraints.HORIZONTAL;
        runPidPanel_.add(runInfoButton_, gridBagConstraints4);
        
        endNowButton_.setToolTipText("Stop job at end of current time step (needs runTimeModifiable)");
        endNowButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        endNowButton_.setText("endNow");
        endNowButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                onEndNow(evt);
            }
        });
        
        gridBagConstraints4 = new java.awt.GridBagConstraints();
        gridBagConstraints4.gridx = 1;
        gridBagConstraints4.gridy = 0;
        gridBagConstraints4.fill = java.awt.GridBagConstraints.HORIZONTAL;
        runPidPanel_.add(endNowButton_, gridBagConstraints4);
        
        endButton_.setToolTipText("Stop job at next dump time step (needs runTimeModifiable)");
        endButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        endButton_.setText("end");
        endButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnEnd(evt);
            }
        });
        
        gridBagConstraints4 = new java.awt.GridBagConstraints();
        gridBagConstraints4.gridx = 2;
        gridBagConstraints4.gridy = 0;
        gridBagConstraints4.fill = java.awt.GridBagConstraints.HORIZONTAL;
        runPidPanel_.add(endButton_, gridBagConstraints4);
        
        killButton_.setToolTipText("Kill running job");
        killButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        killButton_.setText("kill");
        killButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnKill(evt);
            }
        });
        
        gridBagConstraints4 = new java.awt.GridBagConstraints();
        gridBagConstraints4.gridx = 3;
        gridBagConstraints4.gridy = 0;
        gridBagConstraints4.fill = java.awt.GridBagConstraints.HORIZONTAL;
        runPidPanel_.add(killButton_, gridBagConstraints4);
        
        suspButton_.setToolTipText("Suspend running job");
        suspButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        suspButton_.setText("susp");
        suspButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnSusp(evt);
            }
        });
        
        gridBagConstraints4 = new java.awt.GridBagConstraints();
        gridBagConstraints4.gridx = 4;
        gridBagConstraints4.gridy = 0;
        gridBagConstraints4.fill = java.awt.GridBagConstraints.HORIZONTAL;
        runPidPanel_.add(suspButton_, gridBagConstraints4);
        
        contButton_.setToolTipText("Continue suspended job");
        contButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        contButton_.setText("cont");
        contButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnCont(evt);
            }
        });
        
        gridBagConstraints4 = new java.awt.GridBagConstraints();
        gridBagConstraints4.gridx = 5;
        gridBagConstraints4.gridy = 0;
        gridBagConstraints4.fill = java.awt.GridBagConstraints.HORIZONTAL;
        runPidPanel_.add(contButton_, gridBagConstraints4);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.EAST;
        runPanel_.add(runPidPanel_, gridBagConstraints2);
        
        runTableScrollPane_.setVerticalScrollBarPolicy(javax.swing.JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        runTableScrollPane_.setPreferredSize(new java.awt.Dimension(700, 200));
        runTableScrollPane_.setMinimumSize(new java.awt.Dimension(700, 200));
        runTableScrollPane_.setFont(new java.awt.Font("Dialog", 0, 10));
        runTableScrollPane_.setAutoscrolls(true);
        runningTable_.setModel(runningSorter_);
        runningTable_.setFont(new java.awt.Font("Dialog", 0, 10));
        runningTable_.setPreferredSize(new java.awt.Dimension(700, 200));
        runningTable_.setMaximumSize(new java.awt.Dimension(30000, 30000));
        runningTable_.setPreferredScrollableViewportSize(new java.awt.Dimension(0, 0));
        runningTable_.setMinimumSize(new java.awt.Dimension(700, 200));
        runTableScrollPane_.setViewportView(runningTable_);
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 2;
        gridBagConstraints2.gridwidth = 2;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.NORTHWEST;
        gridBagConstraints2.weightx = 1.0;
        gridBagConstraints2.weighty = 1.0;
        runPanel_.add(runTableScrollPane_, gridBagConstraints2);
        
        tabbedPanel_.addTab("runningJobs", runPanel_);
        
        finPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints5;
        
        finPanel_.setBorder(new javax.swing.border.EmptyBorder(new java.awt.Insets(10, 10, 10, 10)));
        finPanel_.setToolTipText("Finished Jobs");
        finPanel_.setPreferredSize(new java.awt.Dimension(730, 161));
        finPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        finAllPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints6;
        
        finAllPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        finAllPanel_.setName("All Processes");
        finAllPanel_.setBackground(new java.awt.Color(204, 204, 204));
        finAllPanel_.setForeground(java.awt.Color.black);
        purgeOldButton_.setToolTipText("Remove old entries");
        purgeOldButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        purgeOldButton_.setText("purge");
        purgeOldButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnPurgeOld(evt);
            }
        });
        
        gridBagConstraints6 = new java.awt.GridBagConstraints();
        gridBagConstraints6.gridx = 0;
        gridBagConstraints6.gridy = 0;
        gridBagConstraints6.fill = java.awt.GridBagConstraints.HORIZONTAL;
        finAllPanel_.add(purgeOldButton_, gridBagConstraints6);
        
        gridBagConstraints5 = new java.awt.GridBagConstraints();
        gridBagConstraints5.gridx = 0;
        gridBagConstraints5.gridy = 0;
        gridBagConstraints5.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints5.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints5.anchor = java.awt.GridBagConstraints.WEST;
        finPanel_.add(finAllPanel_, gridBagConstraints5);
        
        finPidPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints7;
        
        finPidPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        finPidPanel_.setName("Single Procses");
        finPidPanel_.setBackground(new java.awt.Color(204, 204, 204));
        finPidPanel_.setForeground(java.awt.Color.black);
        finInfoButton_.setToolTipText("Display job info");
        finInfoButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        finInfoButton_.setText("Info");
        finInfoButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                onFinInfo(evt);
            }
        });
        
        gridBagConstraints7 = new java.awt.GridBagConstraints();
        gridBagConstraints7.gridx = 0;
        gridBagConstraints7.gridy = 0;
        gridBagConstraints7.fill = java.awt.GridBagConstraints.HORIZONTAL;
        finPidPanel_.add(finInfoButton_, gridBagConstraints7);
        
        removeButton_.setToolTipText("Remove entry");
        removeButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        removeButton_.setText("remove");
        removeButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                onRemove(evt);
            }
        });
        
        gridBagConstraints7 = new java.awt.GridBagConstraints();
        gridBagConstraints7.gridx = 1;
        gridBagConstraints7.gridy = 0;
        gridBagConstraints7.fill = java.awt.GridBagConstraints.HORIZONTAL;
        finPidPanel_.add(removeButton_, gridBagConstraints7);
        
        gridBagConstraints5 = new java.awt.GridBagConstraints();
        gridBagConstraints5.gridx = 1;
        gridBagConstraints5.gridy = 0;
        gridBagConstraints5.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints5.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints5.anchor = java.awt.GridBagConstraints.EAST;
        finPanel_.add(finPidPanel_, gridBagConstraints5);
        
        finTableScrollPane_.setVerticalScrollBarPolicy(javax.swing.JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        finTableScrollPane_.setPreferredSize(new java.awt.Dimension(700, 200));
        finTableScrollPane_.setMinimumSize(new java.awt.Dimension(700, 200));
        finTableScrollPane_.setFont(new java.awt.Font("Dialog", 0, 10));
        finTableScrollPane_.setAutoscrolls(true);
        finishedTable_.setModel(finishedSorter_);
        finishedTable_.setFont(new java.awt.Font("Dialog", 0, 10));
        finishedTable_.setPreferredSize(new java.awt.Dimension(700, 200));
        finishedTable_.setMaximumSize(new java.awt.Dimension(30000, 30000));
        finishedTable_.setPreferredScrollableViewportSize(new java.awt.Dimension(0, 0));
        finishedTable_.setMinimumSize(new java.awt.Dimension(700, 200));
        finTableScrollPane_.setViewportView(finishedTable_);
        
        gridBagConstraints5 = new java.awt.GridBagConstraints();
        gridBagConstraints5.gridx = 0;
        gridBagConstraints5.gridy = 1;
        gridBagConstraints5.gridwidth = 2;
        gridBagConstraints5.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints5.insets = new java.awt.Insets(5, 5, 5, 5);
        gridBagConstraints5.anchor = java.awt.GridBagConstraints.NORTHWEST;
        gridBagConstraints5.weightx = 1.0;
        gridBagConstraints5.weighty = 1.0;
        finPanel_.add(finTableScrollPane_, gridBagConstraints5);
        
        tabbedPanel_.addTab("finishedJobs", finPanel_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.gridwidth = 2;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        getContentPane().add(tabbedPanel_, gridBagConstraints1);
        
        closeButtonPanel_.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.RIGHT));
        
        closeButtonPanel_.setBorder(new javax.swing.border.EtchedBorder());
        closeButtonPanel_.setFont(new java.awt.Font("Dialog", 0, 10));
        myRunsBox_.setToolTipText("Show my jobs only");
        myRunsBox_.setFont(new java.awt.Font("Dialog", 0, 10));
        myRunsBox_.setText("My Jobs");
        myRunsBox_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                myRunsBox_ActionPerformed(evt);
            }
        });
        
        closeButtonPanel_.add(myRunsBox_);
        
        compactBox_.setToolTipText("Compact Representation");
        compactBox_.setFont(new java.awt.Font("Dialog", 0, 10));
        compactBox_.setText("Compact");
        compactBox_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                compactBox_ActionPerformed(evt);
            }
        });
        
        closeButtonPanel_.add(compactBox_);
        
        closeButton_.setToolTipText("Close Panel");
        closeButton_.setFont(new java.awt.Font("Dialog", 0, 10));
        closeButton_.setText("Close");
        closeButton_.addActionListener(new java.awt.event.ActionListener()
        {
            public void actionPerformed(java.awt.event.ActionEvent evt)
            {
                OnClose(evt);
            }
        });
        
        closeButtonPanel_.add(closeButton_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 2;
        gridBagConstraints1.gridwidth = 2;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        getContentPane().add(closeButtonPanel_, gridBagConstraints1);
        
        pack();
        java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        setSize(new java.awt.Dimension(750, 380));
        setLocation((screenSize.width-750)/2,(screenSize.height-380)/2);
    }//GEN-END:initComponents

    private void onRunInfo(java.awt.event.ActionEvent evt)//GEN-FIRST:event_onRunInfo
    {//GEN-HEADEREND:event_onRunInfo
        // Add your handling code here:
        int startRow = runningTable_.getSelectionModel().getMinSelectionIndex();
        int endRow = runningTable_.getSelectionModel().getMaxSelectionIndex();
        if (startRow == endRow)
        {
            displayJobInfo(runningModel_, startRow);
        }
    }//GEN-LAST:event_onRunInfo

    private void onFinInfo(java.awt.event.ActionEvent evt)//GEN-FIRST:event_onFinInfo
    {//GEN-HEADEREND:event_onFinInfo
        // Add your handling code here:
        int startRow = finishedTable_.getSelectionModel().getMinSelectionIndex();
        int endRow = finishedTable_.getSelectionModel().getMaxSelectionIndex();
        if (startRow == endRow)
        {
            displayJobInfo(finishedModel_, startRow);
        }
    }//GEN-LAST:event_onFinInfo

    private void onEndNow(java.awt.event.ActionEvent evt)//GEN-FIRST:event_onEndNow
    {//GEN-HEADEREND:event_onEndNow
        // Add your handling code here:

        int startRow = runningTable_.getSelectionModel().getMinSelectionIndex();
        int endRow = runningTable_.getSelectionModel().getMaxSelectionIndex();

        try
        {
            for (int row = startRow; row <= endRow; row++)
            {
                JobID jobID = runningModel_.mkJobID(row);

                String caseRoot = runningModel_.caseRoot(row);
                String caseName = runningModel_.caseName(row);

                if ((caseRoot.length() == 0) || (caseName.length() == 0))
                {
                    JOptionPane.showConfirmDialog
                    (
                        this,
                        "Process "
                      + runningModel_.toString(jobID)
                      + "\n( "
                      + runningModel_.mkCodeString(row)
                      + " )\n"
                      + "has empty root or case and can not be 'ended'",
                        "Information Dialog",
                        JOptionPane.DEFAULT_OPTION,
                        JOptionPane.INFORMATION_MESSAGE
                    );

                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                }
                else if
                (
                    askConfirm
                    (
                        "Ok to end process "
                      + runningModel_.toString(jobID)
                      + "\n( "
                      + runningModel_.mkCodeString(row)
                      + " )"
                      + "\n(This will stop the job after the next time step)"
                    )
                )
                {
                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                    caseBrowser_.end
                    (
                        jobID,
                        caseRoot,
                        caseName,
                        true                        //  end now
                    );
                    runningModel_.refresh
                    (
                        myRunsBox_.isSelected(),
                        compactBox_.isSelected()
                    );
                    setLicLabel();
                }
                else
                {
                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                }
            }
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch (FoamXSYSError sysErr)
        {
            App.handleAllExceptions(sysErr);
        }

        initTableSize(runningTable_, runningModel_);
        runningTable_.repaint();
        endNowButton_.setEnabled(false);

    }//GEN-LAST:event_onEndNow

    private void compactBox_ActionPerformed(java.awt.event.ActionEvent evt)//GEN-FIRST:event_compactBox_ActionPerformed
    {//GEN-HEADEREND:event_compactBox_ActionPerformed
        // Add your handling code here:
        runningModel_.refresh
        (
            myRunsBox_.isSelected(),
            compactBox_.isSelected()
        );
        setLicLabel();
        initTableSize(runningTable_, runningModel_);
        runningTable_.repaint();

        finishedModel_.refresh
        (
            myRunsBox_.isSelected(),
            compactBox_.isSelected()
        );
        initTableSize(finishedTable_, finishedModel_);
        finishedTable_.repaint();
    }//GEN-LAST:event_compactBox_ActionPerformed

    private void myRunsBox_ActionPerformed(java.awt.event.ActionEvent evt)//GEN-FIRST:event_myRunsBox_ActionPerformed
    {//GEN-HEADEREND:event_myRunsBox_ActionPerformed
        // Add your handling code here:
        runningModel_.refresh
        (
            myRunsBox_.isSelected(),
            compactBox_.isSelected()
        );
        setLicLabel();
        initTableSize(runningTable_, runningModel_);
        runningTable_.repaint();

        finishedModel_.refresh
        (
            myRunsBox_.isSelected(),
            compactBox_.isSelected()
        );
        initTableSize(finishedTable_, finishedModel_);
        finishedTable_.repaint();
    }//GEN-LAST:event_myRunsBox_ActionPerformed

    private void onRemove(java.awt.event.ActionEvent evt)//GEN-FIRST:event_onRemove
    {//GEN-HEADEREND:event_onRemove
        // Add your handling code here:
        try
        {
            int startRow = finishedTable_.getSelectionModel().getMinSelectionIndex();
            int endRow = finishedTable_.getSelectionModel().getMaxSelectionIndex();

            for (int row = startRow; row <= endRow; row++)
            {
                JobID jobID = finishedModel_.mkJobID(row);

                caseBrowser_.purgeFinishedJob(jobID);
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

        finishedTable_.getSelectionModel().clearSelection();
        
        finishedModel_.refresh
        (
            myRunsBox_.isSelected(),
            compactBox_.isSelected()
        );
        initTableSize(finishedTable_, finishedModel_);
        finishedTable_.repaint();

        removeButton_.setEnabled(false);

    }//GEN-LAST:event_onRemove

    private void OnPurgeOld(java.awt.event.ActionEvent evt)//GEN-FIRST:event_OnPurgeOld
    {//GEN-HEADEREND:event_OnPurgeOld
        // Add your handling code here:
        try
        {
            caseBrowser_.purgeFinishedJobs(7);

            finishedModel_.refresh
            (
                myRunsBox_.isSelected(),
                compactBox_.isSelected()
            );

            initTableSize(finishedTable_, finishedModel_);
            finishedTable_.repaint();
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch (FoamXIOError ioErr)
        {
            App.handleAllExceptions(ioErr);
        }
    }//GEN-LAST:event_OnPurgeOld

    private void OnCont(java.awt.event.ActionEvent evt)//GEN-FIRST:event_OnCont
    {//GEN-HEADEREND:event_OnCont
        // Add your handling code here:

        int startRow = runningTable_.getSelectionModel().getMinSelectionIndex();
        int endRow = runningTable_.getSelectionModel().getMaxSelectionIndex();

        try
        {
            for (int row = startRow; row <= endRow; row++)
            {
                JobID jobID = runningModel_.mkJobID(row);

                if
                (
                    askConfirm
                    (
                        "Ok to continue process "
                      + ProcessModel.toString(jobID)
                      + "\n( "
                      + runningModel_.mkCodeString(row)
                      + " )"
                    )
                )
                {
                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                    caseBrowser_.cont(jobID);
                    runningModel_.refresh
                    (
                        myRunsBox_.isSelected(),
                        compactBox_.isSelected()
                    );
                    setLicLabel();
                }
                else
                {
                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                }

            }
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch (FoamXSYSError sysErr)
        {
            App.handleAllExceptions(sysErr);
        }

        initTableSize(runningTable_, runningModel_);
        runningTable_.repaint();

        contButton_.setEnabled(false);

    }//GEN-LAST:event_OnCont

    private void OnSusp(java.awt.event.ActionEvent evt)//GEN-FIRST:event_OnSusp
    {//GEN-HEADEREND:event_OnSusp
        // Add your handling code here:

        int startRow = runningTable_.getSelectionModel().getMinSelectionIndex();
        int endRow = runningTable_.getSelectionModel().getMaxSelectionIndex();

        try
        {
            for (int row = startRow; row <= endRow; row++)
            {
                JobID jobID = runningModel_.mkJobID(row);

                if
                (
                    askConfirm
                    (
                        "Ok to suspend process "
                      + ProcessModel.toString(jobID)
                      + "\n( "
                      + runningModel_.mkCodeString(row)
                      + " )"
                    )
                )
                {
                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                    caseBrowser_.suspend(jobID);
                    runningModel_.refresh
                    (
                        myRunsBox_.isSelected(),
                        compactBox_.isSelected()
                    );
                    setLicLabel();
                }
                else
                {
                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                }
            }
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch (FoamXSYSError sysErr)
        {
            App.handleAllExceptions(sysErr);
        }

        initTableSize(runningTable_, runningModel_);
        runningTable_.repaint();

        suspButton_.setEnabled(false);
        
    }//GEN-LAST:event_OnSusp

    private void OnKill(java.awt.event.ActionEvent evt)//GEN-FIRST:event_OnKill
    {//GEN-HEADEREND:event_OnKill
        // Add your handling code here:

        int startRow = runningTable_.getSelectionModel().getMinSelectionIndex();
        int endRow = runningTable_.getSelectionModel().getMaxSelectionIndex();

        try
        {
            for (int row = startRow; row <= endRow; row++)
            {
                JobID jobID = runningModel_.mkJobID(row);

                if
                (
                    askConfirm
                    (
                        "Ok to kill process "
                      + ProcessModel.toString(jobID)
                      + "\n( "
                      + runningModel_.mkCodeString(row)
                      + " )"
                    )
                )
                {
                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                    caseBrowser_.kill(jobID);
                    runningModel_.refresh
                    (
                        myRunsBox_.isSelected(),
                        compactBox_.isSelected()
                    );
                    setLicLabel();
                    finishedModel_.refresh
                    (
                        myRunsBox_.isSelected(),
                        compactBox_.isSelected()
                    );
                }
                else
                {
                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                }
            }
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch (FoamXSYSError sysErr)
        {
            App.handleAllExceptions(sysErr);
        }

        initTableSize(runningTable_, runningModel_);
        initTableSize(finishedTable_, finishedModel_);
        runningTable_.repaint();

        killButton_.setEnabled(false);
        
    }//GEN-LAST:event_OnKill

    private void OnEnd(java.awt.event.ActionEvent evt)//GEN-FIRST:event_OnEnd
    {//GEN-HEADEREND:event_OnEnd
        // Add your handling code here:

        int startRow = runningTable_.getSelectionModel().getMinSelectionIndex();
        int endRow = runningTable_.getSelectionModel().getMaxSelectionIndex();

        try
        {
            for (int row = startRow; row <= endRow; row++)
            {
                JobID jobID = runningModel_.mkJobID(row);

                String caseRoot = runningModel_.caseRoot(row);
                String caseName = runningModel_.caseName(row);

                if ((caseRoot.length() == 0) || (caseName.length() == 0))
                {
                    JOptionPane.showConfirmDialog
                    (
                        this,
                        "Process "
                      + runningModel_.toString(jobID)
                      + "\n( "
                      + runningModel_.mkCodeString(row)
                      + " )\n"
                      + "has empty root or case and can not be 'ended'",
                        "Information Dialog",
                        JOptionPane.DEFAULT_OPTION,
                        JOptionPane.INFORMATION_MESSAGE
                    );

                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                }
                else if
                (
                    askConfirm
                    (
                        "Ok to end process "
                      + runningModel_.toString(jobID)
                      + "\n( "
                      + runningModel_.mkCodeString(row)
                      + " )"
                      + "\n(This will stop the job when next writing results)"
                    )
                )
                {
                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                    caseBrowser_.end
                    (
                        jobID,
                        runningModel_.caseRoot(row),
                        runningModel_.caseName(row),
                        false                        //  end at next write
                    );
                    runningModel_.refresh
                    (
                        myRunsBox_.isSelected(),
                        compactBox_.isSelected()
                    );
                    setLicLabel();
                }
                else
                {
                    runningTable_.getSelectionModel().removeSelectionInterval
                    (
                        row,
                        row
                    );
                }
            }
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch (FoamXSYSError sysErr)
        {
            App.handleAllExceptions(sysErr);
        }

        initTableSize(runningTable_, runningModel_);
        runningTable_.repaint();

        endButton_.setEnabled(false);
        
    }//GEN-LAST:event_OnEnd

    private void OnPurge(java.awt.event.ActionEvent evt)//GEN-FIRST:event_OnPurge
    {//GEN-HEADEREND:event_OnPurge
        // Add your handling code here:
        try
        {
            caseBrowser_.purgeRunningJobs();
            runningModel_.refresh
            (
                myRunsBox_.isSelected(),
                compactBox_.isSelected()
            );
            setLicLabel();
            initTableSize(runningTable_, runningModel_);
            runningTable_.repaint();
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch (FoamXIOError ioErr)
        {
            App.handleAllExceptions(ioErr);
        }
   }//GEN-LAST:event_OnPurge

    private void OnStatus(java.awt.event.ActionEvent evt)//GEN-FIRST:event_OnStatus
    {//GEN-HEADEREND:event_OnStatus
        // Add your handling code here:
        if
        (
            askConfirm
            (
                "Ok to update status?\n" +
                "(This will contact licensed hosts\n" +
                " and might take some time)"
            )
        )
        {
            // Update status on caseBrowser
            try
            {
                caseBrowser_.refreshJobsLists();
                caseBrowser_.checkRunningJobs();
            }
            catch (FoamXError fxErr)
            {
                App.handleAllExceptions(fxErr);
            }
            catch (FoamXIOError ioErr)
            {
                App.handleAllExceptions(ioErr);
            }
            catch (FoamXSYSError sysErr)
            {
                App.handleAllExceptions(sysErr);
            }

            // Reload data
            runningModel_.refresh
            (
                myRunsBox_.isSelected(),
                compactBox_.isSelected()
            );
            setLicLabel();
            finishedModel_.refresh
            (
                myRunsBox_.isSelected(),
                compactBox_.isSelected()
            );

            runningTable_.updateUI();

            initTableSize(runningTable_, runningModel_);
            runningTable_.repaint();
        }

    }//GEN-LAST:event_OnStatus

    private void OnRead(java.awt.event.ActionEvent evt)//GEN-FIRST:event_OnRead
    {//GEN-HEADEREND:event_OnRead
        // Add your handling code here:

        try
        {
            caseBrowser_.refreshJobsLists();

            runningModel_.refresh
            (
                myRunsBox_.isSelected(),
                compactBox_.isSelected()
            );
            setLicLabel();
            finishedModel_.refresh
            (
                myRunsBox_.isSelected(),
                compactBox_.isSelected()
            );

            runningTable_.updateUI();
            finishedTable_.updateUI();
        }
        catch (FoamXError fxErr)
        {
            App.handleAllExceptions(fxErr);
        }
        catch (FoamXIOError ioErr)
        {
            App.handleAllExceptions(ioErr);
        }

        initTableSize(runningTable_, runningModel_);
        runningTable_.repaint();

        initTableSize(finishedTable_, finishedModel_);
    }//GEN-LAST:event_OnRead



    //--------------------------------------------------------------------------

  private void shapeTypeComboActionPerformed (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_shapeTypeComboActionPerformed
  }//GEN-LAST:event_shapeTypeComboActionPerformed

    //--------------------------------------------------------------------------

  private void OnClose (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnClose

        // Close the dialog.
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
  private javax.swing.JPanel headerPanel_;
  private javax.swing.JLabel licenceLabel_;
  private javax.swing.JTabbedPane tabbedPanel_;
  private javax.swing.JPanel runPanel_;
  private java.awt.Panel runAllPanel_;
  private javax.swing.JButton readButton_;
  private javax.swing.JButton statusButton_;
  private javax.swing.JButton purgeButton_;
  private java.awt.Panel runPidPanel_;
  private javax.swing.JButton runInfoButton_;
  private javax.swing.JButton endNowButton_;
  private javax.swing.JButton endButton_;
  private javax.swing.JButton killButton_;
  private javax.swing.JButton suspButton_;
  private javax.swing.JButton contButton_;
  private javax.swing.JScrollPane runTableScrollPane_;
  private javax.swing.JTable runningTable_;
  private javax.swing.JPanel finPanel_;
  private java.awt.Panel finAllPanel_;
  private javax.swing.JButton purgeOldButton_;
  private java.awt.Panel finPidPanel_;
  private javax.swing.JButton finInfoButton_;
  private javax.swing.JButton removeButton_;
  private javax.swing.JScrollPane finTableScrollPane_;
  private javax.swing.JTable finishedTable_;
  private javax.swing.JPanel closeButtonPanel_;
  private javax.swing.JCheckBox myRunsBox_;
  private javax.swing.JCheckBox compactBox_;
  private javax.swing.JButton closeButton_;
  // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    private class RunningProcessModel
        extends ProcessModel
    {
        public RunningProcessModel(ICaseBrowser caseBrowser)
        {
            super(caseBrowser);
        }

        public void refresh(boolean myJobsOnly, boolean compact)
        {
            try
            {
                ICaseBrowser caseBrowser = getCaseBrowser();
                StringHolder holder = new StringHolder();
                caseBrowser.getUserName(holder);
                
                JobDescriptor[] jobs = caseBrowser.runningJobs();

                setData(myJobsOnly, compact, holder.value, jobs);
                fireTableChanged(new TableModelEvent(this));
            }
            catch (FoamXError fxErr)
            {
                App.handleAllExceptions(fxErr);
            }
        }
    }

    private class FinishedProcessModel
        extends ProcessModel
    {
        public FinishedProcessModel(ICaseBrowser caseBrowser)
        {
            super(caseBrowser);
        }

        public void refresh(boolean myJobsOnly, boolean compact)
        {
            try
            {
                ICaseBrowser caseBrowser = getCaseBrowser();
                StringHolder holder = new StringHolder();
                caseBrowser.getUserName(holder);

                JobDescriptor[] jobs = getCaseBrowser().finishedJobs();

                setData(myJobsOnly, compact, holder.value, jobs);
                fireTableChanged(new TableModelEvent(this));
            }
            catch (FoamXError fxErr)
            {
                App.handleAllExceptions(fxErr);
            }
        }
    }

    //--------------------------------------------------------------------------
}
