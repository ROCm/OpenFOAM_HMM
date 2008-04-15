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

import javax.swing.Icon;

import FoamXServer.CaseBrowser.ICaseBrowser;
import FoamXServer.CaseDescriptor;

import FoamX.App;

public class ContextInfo
    implements FoamX.Util.FoamXTreeRenderer.FoamXTreeItem
{
    //--------------------------------------------------------------------------
    // Symbolic constants for icon indexes.
    protected static final int ROOTCONTEXT = 0;
    protected static final int HOST        = 1;
    protected static final int HOSTOFFLINE = 2;
    protected static final int CASEROOT    = 3;
    protected static final int CASENAME    = 4;
    protected static final int LOCKEDCASE  = 5;
    protected static final int MANAGEDCASE = 6;
    protected static final int ERRORCASE   = 7;
    protected static final int INVALID_REF = 8;

    protected int            type_           = ROOTCONTEXT;
    protected String         name_           = "";
    protected ICaseBrowser   caseBrowser_    = null;
    protected CaseDescriptor caseDescriptor_ = null;
    protected String         caseRoot_       = "";
    protected String         caseName_       = "";
    protected String         app_       = "";
    protected String         statusText_     = "";
    protected boolean        caseLocked_     = false;
    protected boolean        caseManaged_    = false;
    protected boolean        caseError_      = false;

    protected static Icon[] icons_;

    //--------------------------------------------------------------------------

    static
    {
        try
        {
            // Load icons.
            icons_ = new Icon[9];
            icons_[ROOTCONTEXT] =
                App.getResources().getIcon("CaseBrowser.RootContextImage");
            icons_[HOST]        =
                App.getResources().getIcon("CaseBrowser.HostImage");
            icons_[HOSTOFFLINE] =
                App.getResources().getIcon("CaseBrowser.HostOfflineImage");
            icons_[CASEROOT]    =
                App.getResources().getIcon("CaseBrowser.CaseRootImage");
            icons_[CASENAME]    =
                App.getResources().getIcon("CaseBrowser.CaseNameImage");
            icons_[LOCKEDCASE]  =
                App.getResources().getIcon("CaseBrowser.LockedCaseImage");
            icons_[MANAGEDCASE] =
                App.getResources().getIcon("CaseBrowser.ManagedCaseImage");
            icons_[ERRORCASE]   =
                App.getResources().getIcon("CaseBrowser.InvalidRefImage");
            icons_[INVALID_REF] =
                App.getResources().getIcon("CaseBrowser.InvalidRefImage");
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** ContextInfo constructor. */
    public ContextInfo(String name, int type)
    {
        name_ = name;
        type_ = type;
    }

    public ContextInfo
    (
         ICaseBrowser caseBrowser,
         String caseRoot,
         String displayCaseRoot
    )
    {
        name_            = displayCaseRoot;
        type_            = CASEROOT;
        caseBrowser_     = caseBrowser;
        caseRoot_        = caseRoot;
    }

    public ContextInfo(ICaseBrowser caseBrowser, CaseDescriptor caseDescriptor)
    {
        caseBrowser_    = caseBrowser;
        caseDescriptor_ = caseDescriptor;

        name_           = caseDescriptor_.caseName;
        type_           = CASENAME;
        caseRoot_       = caseDescriptor_.rootDir;
        caseName_       = caseDescriptor_.caseName;
        app_       = caseDescriptor_.app;
        caseLocked_     = caseDescriptor_.locked;
        caseManaged_    = caseDescriptor_.managed;
        caseError_      = caseDescriptor_.error;
    }

    //--------------------------------------------------------------------------
    // FoamXTreeItem interface functions.

    public Icon getIcon()
    {
        // See if this case is in error, locked or unmanaged.
        if (type_ == CASENAME && caseError_)
        {
            return icons_[ERRORCASE];
        }
        else if (type_ == CASENAME && caseLocked_)
        {
            return icons_[LOCKEDCASE];
        }
        else if (type_ == CASENAME && caseManaged_)
        {
            return icons_[MANAGEDCASE];
        }
        return icons_[type_];
    }

    public String getText()
    {
        return toString();
    }

    //--------------------------------------------------------------------------

    public String toString()
    {
        return name_;
    }

    public int getType()
    {
        return type_;
    }
    public void setType(int type)
    {
        type_ = type;
    }

    public ICaseBrowser getCaseBrowser()
    {
        return caseBrowser_;
    }
    public void setCaseBrowser(ICaseBrowser caseBrowser)
    {
        caseBrowser_ = caseBrowser;
    }

    public CaseDescriptor getCaseDescriptor()
    {
        return caseDescriptor_;
    }
    public void setCaseDescriptor(CaseDescriptor caseDescriptor)
    {
        caseDescriptor_ = caseDescriptor;
    }

    public String getCaseRoot()
    {
        return caseRoot_;
    }
    public String getCaseName()
    {
        return caseName_;
    }
    public String getAppClass()
    {
        return app_;
    }

    public String getStatusText()
    {
        return statusText_;
    }
    public void setStatusText(String text)
    {
        statusText_ = text;
    }

    public boolean isCaseLocked()
    {
        return caseLocked_;
    }
    public void setCaseLocked(boolean locked)
    {
        caseLocked_ = locked;
    }

    public boolean isCaseError()
    {
        return caseError_;
    }
    public void setCaseError(boolean error)
    {
        caseError_ = error;
    }

    public boolean isCaseManaged()
    {
        return caseManaged_;
    }
    public void setCaseManaged(boolean managed)
    {
        caseManaged_ = managed;
    }

    public boolean isRootContext()
    {
        return type_ == ROOTCONTEXT;
    }
    public boolean isHost()
    {
        return (type_ == HOST || type_ == HOSTOFFLINE);
    }
    public boolean isCaseRoot()
    {
        return type_ == CASEROOT;
    }
    public boolean isCaseName()
    {
        return type_ == CASENAME;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}


