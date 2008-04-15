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

import FoamXServer.CaseServer.ICaseServer;
import FoamXServer.CaseBrowser.ICaseBrowser;

public class CaseSelectionEvent
    extends java.util.EventObject
{
    protected String hostName_;
    protected String caseRoot_;
    protected String caseName_;
    protected ICaseBrowser caseBrowser_;

    /** CaseSelectionEvent constructor. */
    public CaseSelectionEvent
    (
        Object source,
        String hostName,
        String caseRoot,
        String caseName,
        ICaseBrowser caseBrowser
    )
    {
        super(source);
        hostName_    = hostName;
        caseRoot_    = caseRoot;
        caseName_    = caseName;
        caseBrowser_ = caseBrowser;
    }

    public String toString()
    {
        return 
            "CaseSelectionEvent : " + hostName_ + ":" + caseRoot_
            + " " + caseName_;
    }

    public String hostName()
    {
        return hostName_;
    }

    public String caseRoot()
    {
        return caseRoot_;
    }

    public String caseName()
    {
        return caseName_;
    }

    public ICaseBrowser caseBrowser()
    {
        return caseBrowser_;
    }
}



