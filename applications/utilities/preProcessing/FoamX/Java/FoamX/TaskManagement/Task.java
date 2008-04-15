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

public class Task
    extends Object
{
    //--------------------------------------------------------------------------

    protected String taskName_;
    protected int progressValue_;
    protected int progressMin_;
    protected int progressMax_;
    protected TaskStatusBar statusBar_;

    //--------------------------------------------------------------------------
    /** Task constructor. */
    Task(String taskName, int min, int max)
    {
        taskName_      = taskName;
        progressValue_ = min;
        progressMin_   = min;
        progressMax_   = max;
        statusBar_     = null;
    }

    //--------------------------------------------------------------------------

    public void setProgress(int progress)
    {
        progressValue_ = (progress <progressMin_) ? progressMin_ : progress;
        progressValue_ = (progress> progressMax_) ? progressMax_ : progress;

        // Update the status bar.
        if (statusBar_ != null)
        {
            statusBar_.setProgress(progressValue_);
            statusBar_.updateUI();
        }
    }

    //--------------------------------------------------------------------------

    void activate(TaskStatusBar statusBar)
    {
        // Store reference to the status bar.
        statusBar_ = statusBar;

        // Initialise the task status bar with our settings.
        statusBar_.setStatusText(taskName_);
        statusBar_.setProgressMin(progressMin_);
        statusBar_.setProgressMax(progressMax_);
        statusBar_.setProgress(progressValue_);
        //statusBar_.updateUI();
    }

    //--------------------------------------------------------------------------

    void deactivate()
    {
        if (statusBar_ != null)
        {
            // Release reference to status bar.
            statusBar_ = null;
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}


