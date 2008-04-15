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

import javax.swing.table.*;
import javax.swing.event.TableModelListener; 
import javax.swing.event.TableModelEvent; 

import FoamXServer.CaseBrowser.ICaseBrowser;
import FoamXServer.*;

public class ProcessModel
    extends AbstractTableModel
    implements TableModelListener
{
    static public int STAT_INDEX = 0, USER_INDEX = 1, CASE_INDEX = 2;
    static public int HOST_INDEX = 3, PID_INDEX = 4, NPROC_INDEX = 5;
    static public int START_INDEX = 6, END_INDEX = 7, CODE_INDEX = 8;
    static public int ROOT_INDEX = 9, NCNTPROC_INDEX = 10;
    static public int MAX_INDEX = 10;

    static private String[] columnNames_ =
    {
        "Stat", "User", "Case", "Host", "PID", "nCpu",
        "start", "end", "code", "root", "nLic"
    };

    /** Data used to size table columns */
    static private Object[] longValues_ =
    {
        "ENDING",                   // max status width
        "someus",                   // max user width
        "12346789",                 // case
        "12345678",                 // hostname
        "12345678",                 // pid
        "9",                        // nProcs
        "20 Feb ",                  // start data
        "20 Feb ",                  // end data
        "icoFoamJunk",              // code
        "someRoot/aaaa/bbb/ccc",    // root
        "1"                         // nCountedProcs
    };

    private String[] statusNames_;

    protected ICaseBrowser caseBrowser_;

    // Actual jobControl data
    protected Object[][] data_;

    // Sum of all jobs
    protected int nJobs_;
    // Sum of all floatlicensed processes
    protected int nCountedProcs_;


    //--------------------------------------------------------------------------
    /** Creates new empty ProcessModel. Only  for variable initialisation */
    public ProcessModel()
    {
    }


    //--------------------------------------------------------------------------
    /** Creates new ProcessModel */
    public ProcessModel(ICaseBrowser caseBrowser)
    {
        caseBrowser_ = caseBrowser;

        nJobs_ = 0;
        nCountedProcs_ = 0;

        statusNames_ = new String[7];
        statusNames_[FoamXServer.JobStatus.JOB_UNDEFINED.value()] = "OTHR";
        statusNames_[FoamXServer.JobStatus.JOB_LAUNCHING.value()] = "LAUNCH";
        statusNames_[FoamXServer.JobStatus.JOB_RUNNING.value()] = "RUNN";
        statusNames_[FoamXServer.JobStatus.JOB_STOPPING.value()] = "ENDING";
        statusNames_[FoamXServer.JobStatus.JOB_SUSPENDED.value()] = "SUSP";
        statusNames_[FoamXServer.JobStatus.JOB_FINISHED.value()] = "FINI";
        statusNames_[FoamXServer.JobStatus.JOB_ABORTED.value()] = "ABRT";

        refresh(false, false);
    }


    //--------------------------------------------------------------------------
    /** refresh data from caseBrowser. Overloaded */
    public void refresh(boolean myJobsOnly, boolean compact)
    {
        System.out.println("ProcessModel::refresh should never be called");
    }

    public Object[] getLongValues()
    {
        return longValues_;
    }

    public int getNJobs()
    {
        return nJobs_;
    }

    public int getNCountedProcs()
    {
        return nCountedProcs_;
    }


    //--------------------------------------------------------------------------
    /** copy from jobDescriptor[] into data[] array */
    protected void setData
    (
        boolean myJobsOnly,
        boolean compact,
        String userName,
        JobDescriptor[] jobs
    )
    {
        nJobs_ = jobs.length;
        nCountedProcs_ = 0;


        // Count licenced jobs.
        for (int row = 0; row < jobs.length; row++)
        {
            nCountedProcs_ += jobs[row].nCountedProcs;
        }

        // Map of which ones to use.
        int nJobsToUse = 0;
        boolean[] jobsToUse = new boolean[nJobs_];
        for (int i = 0; i < jobs.length; i++)
        {
            jobsToUse[i] = false;
        }

        // Use only user's jobs if specified
        if (myJobsOnly)
        {
            for (int row = 0; row < jobs.length; row++)
            {
                if (jobs[row].userName.equals(userName))
                {
                    jobsToUse[row] = true;
                    nJobsToUse++;
                }
            }
        }
        else
        {
            for (int row = 0; row < jobs.length; row++)
            {
                jobsToUse[row] = true;
                nJobsToUse++;
            }
        }

        // Remove all FoamX stuff if compact
        if (compact)
        {
            for (int row = 0; row < jobs.length; row++)
            {
                if
                (
                    jobsToUse[row]
                 && (
                        (jobs[row].code.equals("FoamXHostBrowser"))
                     || (jobs[row].code.equals("FoamXCaseBrowser"))
                     || (jobs[row].code.equals("FoamXCaseServer"))
                     || (jobs[row].code.equals("FoamXCasePostServer"))
                     || (jobs[row].code.equals("dxFoamExec"))
                     || (jobs[row].code.equals("vxpFoam"))
                     || (jobs[row].code.equals("vxpFoamUser"))
                    )
                )
                {
                    jobsToUse[row] = false;
                    nJobsToUse--;
                }
            }
        }


        // Set data according to jobsToUse

        data_ = new Object[nJobsToUse][MAX_INDEX + 1];
        int i = 0;
        for (int row = 0; row < jobs.length; row++)
        {
            if (jobsToUse[row])
            {
                setDataRow(i++, jobs[row]);
            }
        }
    }

    public ICaseBrowser getCaseBrowser()
    {
        return caseBrowser_;
    }

    private void setDataRow(int row, JobDescriptor job)
    {
        data_[row][STAT_INDEX]     = statusNames_[job.status.value()];
        data_[row][USER_INDEX]     = job.userName;
        data_[row][CASE_INDEX]     = job.caseName;
        data_[row][HOST_INDEX]     = job.jobID.hostName;
        data_[row][PID_INDEX]      = new Integer(job.jobID.processID);
        data_[row][NPROC_INDEX]    = Integer.toString(job.nProcs);
        data_[row][START_INDEX]    = job.startDate;
        data_[row][END_INDEX]      = job.endDate;
        data_[row][CODE_INDEX]     = job.code;
        data_[row][ROOT_INDEX]     = job.rootDir;
        data_[row][NCNTPROC_INDEX] = Integer.toString(job.nCountedProcs);
    }


//--------------------------------------------------------------------------
//
// Implementation of the AbstractTableModel interface
//
//--------------------------------------------------------------------------

    public int getColumnCount()
    {
        return columnNames_.length;
    }

    public int getRowCount()
    {
        return data_.length;
    }

    public String getColumnName(int col)
    {
        return columnNames_[col];
    }

    public Object getValueAt(int row, int col)
    {
        return data_[row][col];
    }


    /*
     * JTable uses this method to determine the default renderer/
     * editor for each cell. (Default would be String)
     */
    public Class getColumnClass(int c)
    {
        if (getRowCount() > 0)
        {
            return getValueAt(0, c).getClass();
        }
        else
        {
            // Empty table. Use String renderer.
            return String.class;
        }
    }

    //
    // Implementation of the TableModelListener interface 
    //

    // By default forward all events to all the listeners. 
    public void tableChanged(TableModelEvent e)
    {
        fireTableChanged(e);
    }

//--------------------------------------------------------------------------

    /** shortcut for accessing root */
    public String caseRoot(int row)
    {
        return (String)getValueAt(row, ROOT_INDEX);
    }

    /** shortcut for accessing casename */
    public String caseName(int row)
    {
        return (String)getValueAt(row, CASE_INDEX);
    }


    /** Create string describing job */
    public String mkCodeString(int row)
    {
        String rootEnding = (String)getValueAt(row, ROOT_INDEX);
        int lastSlash = rootEnding.lastIndexOf('/');
        if (lastSlash == -1)
        {
            lastSlash = 0;
        }

        return
            (String)getValueAt(row, CODE_INDEX) + " ..." +
            rootEnding.substring(lastSlash) + " " +
            (String)getValueAt(row, CASE_INDEX);
    }

    /** get JobID for given entry */
    public JobID mkJobID(int row)
    {
        String hostName =
            (String)getValueAt(row, HOST_INDEX);
        Integer pid2 = (Integer)getValueAt(row, PID_INDEX);

        return new FoamXServer.JobID(hostName, pid2.intValue());
    }

    /** Converted JobID into string */
    public static String toString(JobID jobID)
    {
        return jobID.hostName + Integer.toString(jobID.processID);
    }
}

