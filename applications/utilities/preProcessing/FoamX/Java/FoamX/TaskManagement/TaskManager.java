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
package FoamX.TaskManagement;

import java.util.Stack;
import javax.swing.*;
import java.awt.*;

import FoamX.App;
import FoamX.Exceptions.FoamXException;

public class TaskManager
    extends java.lang.Object
{
    //--------------------------------------------------------------------------

    protected TaskStatusBar statusBar_;
    protected Stack taskStack_;
    protected String defaultText_;

    //--------------------------------------------------------------------------
    /** TaskManager constructor. */
    public TaskManager()
    {
        statusBar_ = new TaskStatusBar();
        taskStack_ = new Stack();

        defaultText_ = "OpenFoamX Version 1.4";

        // Hide the progress bar initially.
        statusBar_.showProgress(false);

        // Set default status text.
        statusBar_.setStatusText(defaultText_);
    }

    //--------------------------------------------------------------------------

    public TaskStatusBar getStatusBar()
    {
        return statusBar_;
    }

    //--------------------------------------------------------------------------

    public Task createTask(String taskName, int progressMin, int progressMax)
    {
        Task newTask = null;

        try
        {
            // Deactivate the previous task if there is one.
            if (taskStack_.size()> 0)
            {
                Task currentTask = (Task)taskStack_.peek();
                currentTask.deactivate();
            }

            // Create a new task object.
            newTask = new Task(taskName, progressMin, progressMax);

            // Add to stack.
            taskStack_.push(newTask);

            // Show busy cursor.
            //FoamXFrame.RootFrame.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

            // Make sure the progress bar is visible.
            //statusBar_.showProgress(true);

            // Activate the new task.
            newTask.activate(statusBar_);
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
            App.handleAllExceptions(ex);
        }

        // Return reference to new task object.
        return newTask;
    }

    //--------------------------------------------------------------------------

    public void endTask()
    {
        try
        {
            if (taskStack_.size() == 0)
            {
                throw new FoamXException("Task count mismatch.");
            }

            // Pop current task off stack.
            Task oldTask = (Task)taskStack_.pop();

            // Deactivate the task.
            oldTask.deactivate();

            // Activate the following task.
            if (taskStack_.size()> 0)
            {
                Task currentTask = (Task)taskStack_.peek();
                currentTask.activate(statusBar_);
            }
            else
            {
                // No more tasks. Hide the progress bar.
                statusBar_.showProgress(false);

                // Restore default status text.
                statusBar_.setStatusText(defaultText_);

                // Restore cursor.
                //FoamXFrame.RootFrame.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            }
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
