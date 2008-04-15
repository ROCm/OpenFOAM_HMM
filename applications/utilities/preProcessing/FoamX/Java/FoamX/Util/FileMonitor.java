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

import java.awt.*;
import java.awt.event.*;
import java.io.*;

import javax.swing.*;
import javax.swing.text.*;
import javax.swing.event.*;

import FoamX.App;
import FoamX.Exceptions.FoamXException;

import FoamXServer.CaseBrowser.ICaseBrowser;


public class FileMonitor implements Runnable
{
    // Objects interested in file status
    protected EventListenerList listenerList_;

    // Casebrowser reference
    protected ICaseBrowser caseBrowser_;

    // Current files to monitor
    protected File[] files_;

    // Delay before checking modification
    protected int sleep_;

    // Last known modification dates
    protected int[] fileDates_;

    // Stop flag
    protected boolean stop_;

    public FileMonitor(ICaseBrowser caseBrowser, File[] files, int sleep)
    {
        // Create listener list.
        listenerList_ = new EventListenerList();

        caseBrowser_ = caseBrowser;

        files_ = files;

        sleep_ = sleep;

        fileDates_ = new int[files.length];

        stop_ = false;

        try
        {
            for(int i = 0; i < files.length; i++)
            {
                String fName = files_[i].getPath();

                fileDates_[i] = caseBrowser_.fileModificationDate(fName);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void run()
    {
        for(;;)
        {
            try
            {
                Thread.currentThread().sleep(sleep_);
            }
            catch (InterruptedException ie)
            {}

            if (stop_)
            {
                return;
            }

            try
            {
                for(int i = 0; i < files_.length; i++)
                {
                    String fName = files_[i].getPath();

                    int currentTime =
                        caseBrowser_.fileModificationDate(fName);

                    if (currentTime != fileDates_[i])
                    {
                        fileDates_[i] = currentTime;

                        fireChangeFile(files_[i]);
                    }
                    else
                    {
                        fileDates_[i] = currentTime;
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

    public void stop()
    {
        stop_ = true;
    }

    //--------------------------------------------------------------------------
    //---- FileListener Interface
    //--------------------------------------------------------------------------

    public void addFileListener(FileListener l)
    {
        listenerList_.add(FileListener.class, l);
    }

    //--------------------------------------------------------------------------

    public void removeFileListener(FileListener l)
    {
        listenerList_.remove(FileListener.class, l);
    }

    //--------------------------------------------------------------------------


    protected void fireOpenFile(File file)
    {
        // Create event object.
        FileEvent evt = new FileEvent(this, file);

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == FileListener.class)
            {
                ((FileListener)listeners[i+1]).fileOpened(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected void fireCloseFile(File file)
    {
        // Create event object.
        FileEvent evt = new FileEvent(this, file);

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == FileListener.class)
            {
                ((FileListener)listeners[i+1]).fileClosed(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected void fireChangeFile(File file)
    {
        // Create event object.
        FileEvent evt = new FileEvent(this, file);

        // Process the listeners last to first, notifying those that are
        // interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == FileListener.class)
            {
                ((FileListener)listeners[i+1]).fileChanged(evt);
            }
        }
    }

    //--------------------------------------------------------------------------
}
