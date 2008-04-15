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

import FoamX.App;
import FoamXServer.ApplicationDescriptor;

public class CreateCaseDialog
    extends javax.swing.JDialog
{
    //--------------------------------------------------------------------------

    protected boolean wasCancelled_;
    protected String caseClass_;
    protected String caseRoot_;
    protected ApplicationDescriptor[] applicationes_;

    //--------------------------------------------------------------------------
    /** CreateCaseDialog constructor. */
    public CreateCaseDialog
    (
        java.awt.Frame parent,
        ApplicationDescriptor[] applicationes,
        String[] caseRoots
    )
    {
        super(parent, true);      // Modal.

        wasCancelled_       = false;
        applicationes_ = applicationes;

        // Initialise the default GUI.
        initComponents();

        // Initialise the application class combo box.
        for (int i = 0; i <applicationes_.length; i++)
        {
            appCombo_.addItem(applicationes_[i].name);
        }

        // Initialise the case roots combo box.
        for (int i = 0; i <caseRoots.length; i++)
        {
            rootCombo_.addItem(caseRoots[i]);
        }

        // Repack everything.
        pack();
    }

    //--------------------------------------------------------------------------

//    /** Initialise the application class combo box. */
//    protected void buildAppsMenu
//    (
//        JComboBox appsMenu,
//        ApplicationDescriptor[] applicationes
//    )
//    {
//        // map from category to menu
//        Hashtable appsMenuMap = new Hashtable(5);
//
//        for (int i = 0; i <applicationes_.length; i++)
//        {
//            JComboBox menu = appsMenu;
//
//            String appName = applicationes_[i].name;
//            String appCat = applicationes_[i].category;
//            if ((appCat != null) && (appCat.length() > 0))
//            {
//                // Has category. Build hierarchical submenus by splitting
//                // category at '/'
//                // Set menu to be submenu.
//                StringTokenizer tok =
//                    new StringTokenizer(appCat, "/");
//
//                String catKey = new String();
//
//                while (tok.hasMoreTokens())
//                {
//                    String category = tok.nextToken();
//                    catKey = catKey + category;
//                    if (appsMenuMap.containsKey(catKey))
//                    {
//                        // Menu already exists for this category
//                        menu = (JComboBox)appsMenuMap.get(catKey);
//                    }
//                    else
//                    {
//                        // No menu yet for this category
//                        JComboBox catMenu = new JComboBox(category);
//                        catMenu.setFont(appsMenu.getFont());
//                        menu.add(catMenu);
//                        appsMenuMap.put(catKey, catMenu);
//                        menu = catMenu;
//                    }
//                }
//            }
//
//            // Handle item itself
//
//            // Add a menu item for this utility.
//            JComboItem comboItem = new JComboItem(appName);
//            comboItem.setFont(appsMenu.getFont());
//
//            // Add to menu (is either toplevel or category menu)
//            menu.add(comboItem);
//        }
//
//    }

    //--------------------------------------------------------------------------

    public boolean wasCancelled()
    {
        return wasCancelled_;
    }

    //--------------------------------------------------------------------------

    public String getCaseClass()
    {
        return caseClass_;
    }
    public void   setCaseClass(String app)
    {
        // Update the application class combo box.
        appCombo_.setSelectedItem(app);
    }

    //--------------------------------------------------------------------------

    public String getCaseRoot()
    {
        return caseRoot_;
    }
    public void   setCaseRoot(String root)
    {
        // Update the root combo box.
        rootCombo_.setSelectedItem(root);
    }

    //--------------------------------------------------------------------------

    public String getCaseName()
    {
        return caseNameTextField_.getText();
    }
    public void   setCaseName(String name)
    {
        caseNameTextField_.setText(name);
    }

    //--------------------------------------------------------------------------
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the FormEditor.
     */
    private void initComponents() {//GEN-BEGIN:initComponents
        mainPanel_ = new javax.swing.JPanel();
        label1_ = new javax.swing.JLabel();
        appCombo_ = new javax.swing.JComboBox();
        label2_ = new javax.swing.JLabel();
        rootCombo_ = new javax.swing.JComboBox();
        label3_ = new javax.swing.JLabel();
        caseNameTextField_ = new javax.swing.JTextField();
        buttonPanel_ = new javax.swing.JPanel();
        buttonOK_ = new javax.swing.JButton();
        buttonCancel_ = new javax.swing.JButton();
        
        getContentPane().setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints1;
        
        setTitle("Create New Case");
        setModal(true);
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                closeDialog(evt);
            }
        });
        
        mainPanel_.setLayout(new java.awt.GridBagLayout());
        java.awt.GridBagConstraints gridBagConstraints2;
        
        mainPanel_.setPreferredSize(new java.awt.Dimension(350, 100));
        label1_.setText("Class");
        label1_.setForeground(java.awt.Color.black);
        label1_.setFont(new java.awt.Font("Dialog", 0, 10));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(10, 10, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        mainPanel_.add(label1_, gridBagConstraints2);
        
        appCombo_.setFont(new java.awt.Font("Dialog", 0, 10));
        appCombo_.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                OnAppClassComboChanged(evt);
            }
        });
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 0;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(10, 5, 5, 10);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        gridBagConstraints2.weightx = 1.0;
        mainPanel_.add(appCombo_, gridBagConstraints2);
        
        label2_.setText("Case Root");
        label2_.setForeground(java.awt.Color.black);
        label2_.setFont(new java.awt.Font("Dialog", 0, 10));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 1;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 10, 5, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        mainPanel_.add(label2_, gridBagConstraints2);
        
        rootCombo_.setFont(new java.awt.Font("Dialog", 0, 10));
        rootCombo_.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                OnRootComboChanged(evt);
            }
        });
        
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 1;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 5, 10);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        gridBagConstraints2.weightx = 1.0;
        mainPanel_.add(rootCombo_, gridBagConstraints2);
        
        label3_.setText("Case Name");
        label3_.setForeground(java.awt.Color.black);
        label3_.setFont(new java.awt.Font("Dialog", 0, 10));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 0;
        gridBagConstraints2.gridy = 2;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 10, 10, 5);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        mainPanel_.add(label3_, gridBagConstraints2);
        
        caseNameTextField_.setFont(new java.awt.Font("Default", 0, 10));
        gridBagConstraints2 = new java.awt.GridBagConstraints();
        gridBagConstraints2.gridx = 1;
        gridBagConstraints2.gridy = 2;
        gridBagConstraints2.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints2.insets = new java.awt.Insets(5, 5, 10, 10);
        gridBagConstraints2.anchor = java.awt.GridBagConstraints.WEST;
        gridBagConstraints2.weightx = 1.0;
        mainPanel_.add(caseNameTextField_, gridBagConstraints2);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 0;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.NORTH;
        gridBagConstraints1.weightx = 1.0;
        gridBagConstraints1.weighty = 1.0;
        getContentPane().add(mainPanel_, gridBagConstraints1);
        
        buttonPanel_.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.RIGHT));
        
        buttonPanel_.setBorder(new javax.swing.border.EtchedBorder());
        buttonOK_.setFont(new java.awt.Font("Dialog", 0, 10));
        buttonOK_.setText("OK");
        buttonOK_.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                OnOK(evt);
            }
        });
        
        buttonPanel_.add(buttonOK_);
        
        buttonCancel_.setFont(new java.awt.Font("Dialog", 0, 10));
        buttonCancel_.setText("Cancel");
        buttonCancel_.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                OnCancel(evt);
            }
        });
        
        buttonPanel_.add(buttonCancel_);
        
        gridBagConstraints1 = new java.awt.GridBagConstraints();
        gridBagConstraints1.gridx = 0;
        gridBagConstraints1.gridy = 1;
        gridBagConstraints1.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints1.anchor = java.awt.GridBagConstraints.SOUTH;
        getContentPane().add(buttonPanel_, gridBagConstraints1);
        
        pack();
        java.awt.Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        setSize(new java.awt.Dimension(400, 300));
        setLocation((screenSize.width-400)/2,(screenSize.height-300)/2);
    }//GEN-END:initComponents

    //--------------------------------------------------------------------------

  private void OnRootComboChanged (java.awt.event.ItemEvent evt)
    {//GEN-FIRST:event_OnRootComboChanged
        // Get the selected case root.
        if (evt.getStateChange() == java.awt.event.ItemEvent.SELECTED)
        {
            caseRoot_ = (String)evt.getItem();
        }
  }//GEN-LAST:event_OnRootComboChanged

    //--------------------------------------------------------------------------

  private void OnAppClassComboChanged (java.awt.event.ItemEvent evt)
    {//GEN-FIRST:event_OnAppClassComboChanged
        // Get the selected application class.
        if (evt.getStateChange() == java.awt.event.ItemEvent.SELECTED)
        {
            caseClass_ = (String)evt.getItem();
        }
  }//GEN-LAST:event_OnAppClassComboChanged

    //--------------------------------------------------------------------------

  private void OnCancel (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnCancel
        wasCancelled_ = true;
        setVisible(false);
        dispose();
  }//GEN-LAST:event_OnCancel

    //--------------------------------------------------------------------------

  private void OnOK (java.awt.event.ActionEvent evt)
    {//GEN-FIRST:event_OnOK

        //JTextField edit = (JTextField)rootCombo_.getEditor().getEditorComponent();
        //edit.setText(caseRoot_);

        // Only close if we have a case name entered.
        if (getCaseName().trim().length()> 0)
        {
            wasCancelled_ = false;
            setVisible(false);
            dispose();
        }
        else
        {
            App.getRootFrame().getToolkit().getDefaultToolkit().beep();
        }

  }//GEN-LAST:event_OnOK

    //--------------------------------------------------------------------------

    /** Closes the dialog */
    private void closeDialog(java.awt.event.WindowEvent evt)
    {//GEN-FIRST:event_closeDialog
        wasCancelled_ = true;
        setVisible(false);
        dispose();
    }//GEN-LAST:event_closeDialog

    //--------------------------------------------------------------------------

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel mainPanel_;
    private javax.swing.JLabel label1_;
    private javax.swing.JComboBox appCombo_;
    private javax.swing.JLabel label2_;
    private javax.swing.JComboBox rootCombo_;
    private javax.swing.JLabel label3_;
    private javax.swing.JTextField caseNameTextField_;
    private javax.swing.JPanel buttonPanel_;
    private javax.swing.JButton buttonOK_;
    private javax.swing.JButton buttonCancel_;
    // End of variables declaration//GEN-END:variables

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

}




