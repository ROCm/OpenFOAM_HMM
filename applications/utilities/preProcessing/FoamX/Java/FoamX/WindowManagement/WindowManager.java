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
package FoamX.WindowManagement;

import java.util.*;
import java.awt.Rectangle;
import java.awt.Point;
import java.awt.event.*;
import javax.swing.*;

import FoamX.App;

public class WindowManager
    implements FoamX.ActionManagement.IActionProvider
{
    //--------------------------------------------------------------------------

    protected JDesktopPane desktopPanel_;      // Reference to the desktop panel object.
    protected JPopupMenu contextMenu_;
    protected Hashtable windowMap_;         // Map of FoamXInternalFrame objects.
    protected int windowCount_;
    protected Point windowPosition_;
    protected int defaultFrameWidth_;
    protected int defaultFrameHeight_;

    //--------------------------------------------------------------------------
    // Window manager actions.
    private Action[] actions_ =
    {
        new ArrangeWindowsCascadeAction(),
        new ArrangeWindowsHorizontalAction(),
        new ArrangeWindowsVerticalAction()
    };

    //--------------------------------------------------------------------------
    /** ModuleWindowManager constructor. */
    public WindowManager(JDesktopPane desktopPanel)
    {
        try
        {
            desktopPanel_   = desktopPanel;
            contextMenu_    = new javax.swing.JPopupMenu();
            windowMap_      = new Hashtable(10);
            windowCount_    = 0;
            windowPosition_ = new Point(0, 0);

            contextMenu_.setFont(new java.awt.Font("Courier", 0, 10));

            // Listen out for right mouse clicks.
            desktopPanel_.addMouseListener
            (
                new java.awt.event.MouseAdapter()
                {
                    public void mouseClicked(java.awt.event.MouseEvent evt)
                    {
                        OnMouseClick(evt);
                    }
                }
            );

            // Get the default frame size.
            defaultFrameWidth_  = 400;
            defaultFrameHeight_ = 350;
            if (App.getOptions().containsKey("WindowManager.InternalFrameWidth"))
            {
                defaultFrameWidth_ = Integer.parseInt(App.getOptions().getProperty("WindowManager.InternalFrameWidth"));
            }
            if (App.getOptions().containsKey("WindowManager.InternalFrameHeight"))
            {
                defaultFrameHeight_ = Integer.parseInt(App.getOptions().getProperty("WindowManager.InternalFrameHeight"));
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Create a new window and adds it to the desktop.
     * @return The id of the new window.
     * @param component A reference to the component to be hosted by the new
     *        window.
     */
    public int createWindow(javax.swing.JComponent component)
    {
        int windowID = -1;

        try
        {
            // Allocate a new Window ID.
            windowID = windowCount_++;

            // Create a new FoamXInternalFrame object and add the component.
            FoamXInternalFrame frame = new FoamXInternalFrame
            (
                windowID,
                component
            );

            // Add new window to the map.
            windowMap_.put(new Integer(windowID), frame);

            // Add new window to the desktop.
            desktopPanel_.add(frame, javax.swing.JLayeredPane.DEFAULT_LAYER);

            // Reset initial position if desktop panel is smaller than required.
            if (windowPosition_.getX()>= 50)   //desktopPanel_.getWidth())
            {
                windowPosition_.setLocation(0, 0);
            }

            int width = defaultFrameWidth_;
            int height = defaultFrameHeight_;

            // Use component provided sizes
            java.awt.Dimension componentSize = component.getPreferredSize();
            if (componentSize != null)
            {
                width = (int)componentSize.getWidth();
                height = (int)componentSize.getHeight();
            }

            // Clip to desktop size
            width = Math.min(width, desktopPanel_.getWidth());
            height = Math.min(height, desktopPanel_.getHeight());

            // Set the frame size and position.
            frame.setBounds
            (
                (int)windowPosition_.getX(),
                (int)windowPosition_.getY(),
                width,
                height
            );
            windowPosition_.translate(10, 10);

            // Subscribe to the FoamXInternalFrame object's frameClosed event.
            frame.addInternalFrameListener
            (
                new javax.swing.event.InternalFrameAdapter()
                {
                    public void internalFrameClosed
                    (
                        javax.swing.event.InternalFrameEvent evt
                    )
                    {
                        FoamXInternalFrame evtFrame =
                            (FoamXInternalFrame)evt.getSource();
                        if (evtFrame != null)
                        {
                            onFrameClosed(evtFrame.getWindowID());
                        }
                    }
                }
            );

            // Finally show the window.
            frame.toFront();
            frame.setVisible(true);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
            windowID = -1;
        }

        return windowID;
    }

    //--------------------------------------------------------------------------

    public void closeWindow(int windowID)
    {
        try
        {
            if (windowMap_.containsKey(new Integer(windowID)))
            {
                FoamXInternalFrame frame = (FoamXInternalFrame)windowMap_.get(new Integer(windowID));

                // Remove window from the desktop.
                desktopPanel_.remove(frame);
                desktopPanel_.updateUI();          // TODO : Check whether we need 

                // Close the window. This will fire a Frame Closed event which
                // we handle in onFrameClosed where we tidy up the window map.
                frame.setClosed(true);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void showWindow(int windowID)
    {
        try
        {
            if (windowMap_.containsKey(new Integer(windowID)))
            {
                FoamXInternalFrame frame = (FoamXInternalFrame)windowMap_.get(new Integer(windowID));

                // See if this window is iconified.
                if (frame.isIcon()) frame.setIcon(false);
                frame.setVisible(true);
                frame.toFront();
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void hideWindow(int windowID)
    {
        try
        {
            if (windowMap_.containsKey(new Integer(windowID)))
            {
                FoamXInternalFrame frame = (FoamXInternalFrame)windowMap_.get(new Integer(windowID));

                // Hide window.
                frame.setVisible(false);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public FoamXInternalFrame getWindowFrame(int windowID)
    {
        FoamXInternalFrame frame = null;

        try
        {
            // Check for valid window ID.
            if (windowMap_.containsKey(new Integer(windowID)))
            {
                frame = (FoamXInternalFrame)windowMap_.get(new Integer(windowID));
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return frame;
    }

    //--------------------------------------------------------------------------

    protected void onFrameClosed(int windowId)
    {
        try
        {
            if (windowMap_.containsKey(new Integer(windowId)))
            {
                // Remove window from map.
                windowMap_.remove(new Integer(windowId));
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    protected void OnMouseClick(java.awt.event.MouseEvent evt)
    {
        try
        {
            // Show context menu for the desktop.
            //if (evt.isPopupTrigger())
            if (evt.getModifiers() == evt.BUTTON3_MASK)
            {
                // Construct context menu.
                contextMenu_.removeAll();

                JMenuItem mnuProperties = new JMenuItem();
                mnuProperties.setFont(new java.awt.Font("Courier", 0, 10));
                mnuProperties.setText("Properties...");
                contextMenu_.add(mnuProperties);

                // Display the context menu.
                contextMenu_.show(desktopPanel_, evt.getX(), evt.getY());
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    protected void arrangeWindowsCascade()
    {
        try
        {
            int x = 0;

            Enumeration enum = windowMap_.elements();
            while (enum.hasMoreElements())
            {
                FoamXInternalFrame frame = (FoamXInternalFrame)enum.nextElement();

                frame.setBounds(x, x, defaultFrameWidth_, defaultFrameHeight_);
                x += 20;
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //---- IActionProvider Interface Methods
    //--------------------------------------------------------------------------

    public String getPrefix()
    {
        return "WindowManager";
    }
    public javax.swing.Action[] getActions()
    {
        return actions_;
    }

    //--------------------------------------------------------------------------
    //---- Action Classes
    //--------------------------------------------------------------------------

    class ArrangeWindowsCascadeAction
    extends AbstractAction
    {
        ArrangeWindowsCascadeAction()
        {
            putValue(Action.NAME, "ArrangeCascade");
            putValue(Action.SMALL_ICON, App.getResources().getIcon("NewPatchPhysicalTypeImage"));
            putValue(Action.SHORT_DESCRIPTION, "Arrange windows cascade");
            putValue(Action.LONG_DESCRIPTION, "Arrange windows cascade");
        }

        public void actionPerformed(ActionEvent e)
        {
            arrangeWindowsCascade();
        }
    }

    //--------------------------------------------------------------------------

    class ArrangeWindowsHorizontalAction
    extends AbstractAction
    {
        ArrangeWindowsHorizontalAction()
        {
            putValue(Action.NAME, "ArrangeHorizontal");
            putValue(Action.SMALL_ICON, App.getResources().getIcon("NewPatchPhysicalTypeImage"));
            putValue(Action.SHORT_DESCRIPTION, "Arrange windows horizontally");
            putValue(Action.LONG_DESCRIPTION, "Arrange windows horizontally");
        }

        public void actionPerformed(ActionEvent e)
        {
            arrangeWindowsCascade();
        }
    }

    //--------------------------------------------------------------------------

    class ArrangeWindowsVerticalAction
    extends AbstractAction
    {
        ArrangeWindowsVerticalAction()
        {
            putValue(Action.NAME, "ArrangeVertical");
            putValue(Action.SMALL_ICON, App.getResources().getIcon("NewPatchPhysicalTypeImage"));
            putValue(Action.SHORT_DESCRIPTION, "Arrange windows vertically");
            putValue(Action.LONG_DESCRIPTION, "Arrange windows vertically");
        }

        public void actionPerformed(ActionEvent e)
        {
            arrangeWindowsCascade();
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}




