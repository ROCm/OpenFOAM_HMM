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

import java.util.*;
import javax.swing.*;
import java.awt.event.*;

import FoamX.App;
import FoamX.Exceptions.FoamXException;
import FoamX.CaseManagement.CaseBrowserPanel;

import FoamXServer.CaseBrowser.ICaseBrowser;

public class CaseChooserDlg
    extends javax.swing.JDialog
    implements CaseSelectionListener
{
    //--------------------------------------------------------------------------

    protected CaseBrowserPanel caseBrowserPanel_;

    protected String hostName_;
    protected String caseRoot_;
    protected String caseName_;
    protected ICaseBrowser caseBrowser_;

    private static final int DEFAULT_WIDTH = 300;
    private static final int DEFAULT_HEIGHT = 400;

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    /** CaseChooserDlg constructor. */
    public CaseChooserDlg(java.awt.Frame parent, int selectMode)
    {
        super(parent, "Select Case", true);

        try
        {

            // Get current caseBrowserModel from existing caseBrowserPanel
            CaseManager caseManager = App.getCaseManager();

            if (caseManager == null)
            {
                throw new FoamXException
                (
                    "No caseManager in top-level App."
                );
            }

            CaseBrowserModel model =
                caseManager.getCaseBrowser().getModel();

            // Open new caseBrowserPanel with existing model in
            // case-selection mode.
            caseBrowserPanel_ =
                new CaseBrowserPanel(model, selectMode);

            initComponents();

            // Set title
            if (selectMode == CaseBrowserPanel.SELECT_ROOT_MODE)
            {
                setTitle("Select Root");
            }
            else if (selectMode == CaseBrowserPanel.SELECT_CASE_MODE)
            {
                setTitle("Select Case");
            }
            else if (selectMode == CaseBrowserPanel.OPEN_CASE_MODE)
            {
                setTitle("Select Case");
            }

            // Clear all vars
            clearSelection();

            // Get warned of any selection
            caseBrowserPanel_.addCaseSelectionListener(this);
        }
        catch (FoamXException ex)
        {
            App.handleAllExceptions(ex);
        }
    }



    /** GUI initialisation */
    private void initComponents()
    {
        getContentPane().add(caseBrowserPanel_);

        addWindowListener
        (
            new java.awt.event.WindowAdapter()
            {
                public void windowClosing(java.awt.event.WindowEvent evt)
                {
                    closeDialog(evt);
                }
            }
        );

        caseBrowserPanel_.setPreferredSize
        (
            new java.awt.Dimension(DEFAULT_WIDTH, DEFAULT_HEIGHT)
        );
        pack();

        java.awt.Dimension screenSize =
            java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        setLocation
        (
            (screenSize.width-getSize().width)/2,
            (screenSize.height-getSize().height)/2
        );
    }

    private void closeDialog(java.awt.event.WindowEvent evt)
    {
        clearSelection();
        exitWindow();
    }

    private void clearSelection()
    {
        hostName_ = null;
        caseRoot_ = null;
        caseName_ = null;
        caseBrowser_ = null;
    }

    private void exitWindow()
    {
        //caseBrowserPanel_.shutdown();
        setVisible(false);
        dispose();
    }

    //--------------------------------------------------------------------------
    //---- public accessors
    //--------------------------------------------------------------------------

    /** get selected host */
    public String getHostName()
    {
        return hostName_;
    }

    /** get selected root */
    public String getCaseRoot()
    {
        return caseRoot_;
    }

    /** get selected case name */
    public String getCaseName()
    {
        return caseName_;
    }

    /** get selected ICaseBrowser */
    public ICaseBrowser getCaseBrowser()
    {
        return caseBrowser_;
    }

    //--------------------------------------------------------------------------
    //---- CaseSelectionListener Interface
    //--------------------------------------------------------------------------


    public void hostSelected(CaseSelectionEvent evt)
    {
        hostName_ = evt.hostName();
        caseRoot_ = evt.caseRoot();
        caseName_ = evt.caseName();
        caseBrowser_ = evt.caseBrowser();
    }

    public void rootSelected(CaseSelectionEvent evt)
    {
        hostName_ = evt.hostName();
        caseRoot_ = evt.caseRoot();
        caseName_ = evt.caseName();
        caseBrowser_ = evt.caseBrowser();
        exitWindow();
    }

    public void caseSelected(CaseSelectionEvent evt)
    {
        hostName_ = evt.hostName();
        caseRoot_ = evt.caseRoot();
        caseName_ = evt.caseName();
        caseBrowser_ = evt.caseBrowser();
        exitWindow();
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}







