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

import FoamX.App;
import FoamX.Exceptions.FoamXException;
import FoamX.CaseManagement.CaseBrowserPanel;

import FoamXServer.CaseBrowser.ICaseBrowser;

public class HostChooserDlg
    extends javax.swing.JDialog
    implements CaseSelectionListener
{
    //--------------------------------------------------------------------------

    protected CaseBrowserPanel caseBrowserPanel_;

    protected String hostName_;

    private static final int DEFAULT_WIDTH = 300;
    private static final int DEFAULT_HEIGHT = 400;

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    /** HostChooserDlg constructor. */
    public HostChooserDlg(java.awt.Frame parent)
    {
        super(parent, "Select Host", true);

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
                new CaseBrowserPanel(model, CaseBrowserPanel.SELECT_HOST_MODE);

            initComponents();

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

    /** get selected root */
    public String getHostName()
    {
        return hostName_;
    }

    //--------------------------------------------------------------------------
    //---- CaseSelectionListener Interface
    //--------------------------------------------------------------------------


    public void hostSelected(CaseSelectionEvent evt)
    {
        hostName_ = evt.hostName();
        exitWindow();
    }

    public void rootSelected(CaseSelectionEvent evt)
    {
    }

    public void caseSelected(CaseSelectionEvent evt)
    {
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}







