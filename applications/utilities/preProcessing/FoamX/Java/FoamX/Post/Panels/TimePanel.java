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
package FoamX.Post.Panels;

import java.awt.Container;
import java.awt.GridLayout;
import java.awt.Dimension;
import java.awt.event.*;
import java.awt.event.ActionEvent;
import java.awt.Graphics;
import java.awt.GraphicsConfiguration;
import java.io.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;

import com.sun.j3d.utils.universe.*;
import javax.media.j3d.*;
import javax.vecmath.*;

import FoamXServer.*;
import FoamXServer.CasePostServer.ICasePostServer;

import FoamX.App;
import FoamX.Post.Panels.*;
import FoamX.Post.PostApp;

public class TimePanel extends javax.swing.JPanel
{
    //--------------------------------------------------------------------------

    protected ICasePostServer casePostServer_;
    protected PostApp pp_;
    protected String caseRoot_;
    protected String caseName_;
    protected EventListenerList listenerList_;


    //--------------------------------------------------------------------------
    /** Creates new form TimePanel */
    public TimePanel
    (
        ICasePostServer casePostServer,
        PostApp pp,
        String caseRoot,
        String caseName
    )
    {
        super();

        // Create listener list.
        listenerList_ = new EventListenerList();

        casePostServer_ = casePostServer;
        pp_ = pp;
        caseRoot_ = caseRoot;
        caseName_ = caseName;

        initComponents();

        updateButton_ActionPerformed(new ActionEvent(this, 0, "dummy"));
    }

    //--------------------------------------------------------------------------

    private void initComponents()
    {
        updateButton_ = new javax.swing.JButton("Read Times");
        timeBox_ = new javax.swing.JComboBox();

        // Globalpanel
        this.setLayout(new GridLayout(1, 2));
        this.add(updateButton_);
        this.add(timeBox_);

        // Listeners
        updateButton_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    updateButton_ActionPerformed(evt);
                }
            }
        );
        timeBox_.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    timeBox_ActionPerformed(evt);
                }
            }
        );

//        java.awt.Dimension screenSize =
//            java.awt.Toolkit.getDefaultToolkit().getScreenSize();
//        //setSize(new java.awt.Dimension(300, 200));
//        setLocation((screenSize.width-480)/2,(screenSize.height-200)/2);
    }//GEN-END:initComponents


    private void updateButton_ActionPerformed(ActionEvent evt)
    {
        String[] timeSteps = casePostServer_.availableTimeSteps();
        timeBox_.removeAllItems();
        for(int i = 0; i < timeSteps.length; i++)
        {
            timeBox_.addItem(timeSteps[i]);
        }
        //timeBox_ActionPerformed(new ActionEvent(this, 0, "dummy"));
    }

    private void timeBox_ActionPerformed(ActionEvent evt)
    {
        try
        {
            if (timeBox_.getSelectedItem() != null)
            {
                String timeName = (String)timeBox_.getSelectedItem();

                int timeIndex = timeBox_.getSelectedIndex();

                casePostServer_.setTime(timeName, timeIndex);

                fireTimeChanged(timeIndex, Double.parseDouble(timeName));
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
    }

    //--------------------------------------------------------------------------
    //---- TimeChangeListener Methods
    //--------------------------------------------------------------------------

    public void addTimeChangeListener(TimeChangeListener l)
    {
        listenerList_.add(TimeChangeListener.class, l);
    }

    //--------------------------------------------------------------------------

    public void removeTimeChangeListener(TimeChangeListener l)
    {
        listenerList_.remove(TimeChangeListener.class, l);
    }

    //--------------------------------------------------------------------------

    protected void fireTimeChanged(int timeIndex, double timeValue)
    {
        // Create event object.
        TimeChangeEvent evt = new TimeChangeEvent(this, timeIndex, timeValue);

        // Process the listeners last to first, notifying those that
        // are interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == TimeChangeListener.class)
            {
                ((TimeChangeListener)listeners[i+1]).timeChanged(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

  // Variables declaration - do not modify//GEN-BEGIN:variables
  private javax.swing.JButton updateButton_;
  private javax.swing.JComboBox timeBox_;

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



