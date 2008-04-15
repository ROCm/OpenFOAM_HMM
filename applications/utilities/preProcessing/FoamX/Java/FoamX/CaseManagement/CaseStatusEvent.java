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
import FoamXServer.CasePostServer.ICasePostServer;
import FoamXServer.CaseBrowser.ICaseBrowser;

public class CaseStatusEvent
    extends java.util.EventObject
{
    protected String caseRoot_;
    protected String caseName_;
    protected ICaseServer caseServer_;
    protected ICasePostServer casePostServer_;
    protected ICaseBrowser caseBrowser_;

    /** CaseStatusEvent constructor. */
    public CaseStatusEvent
    (
        Object source,
        String caseRoot,
        String caseName,
        ICaseServer caseServer,
        ICasePostServer casePostServer,
        ICaseBrowser caseBrowser
    )
    {
        super(source);
        caseRoot_       = caseRoot;
        caseName_       = caseName;
        caseServer_     = caseServer;
        casePostServer_ = casePostServer;
        caseBrowser_    = caseBrowser;
    }

    /** CaseStatusEvent constructor. */
    public CaseStatusEvent
    (
        Object source,
        String caseRoot,
        String caseName,
        ICaseServer caseServer
    )
    {
        this(source, caseRoot, caseName, caseServer, null, null);
    }

    public CaseStatusEvent
    (
        Object source,
        String caseRoot,
        String caseName,
        ICasePostServer casePostServer
    )
    {
        this(source, caseRoot, caseName, null, casePostServer, null);
    }

    public CaseStatusEvent
    (
        Object source,
        String caseRoot,
        String caseName
    )
    {
        this(source, caseRoot, caseName, null, null, null);
    }

    public String toString()
    {
        return "CaseStatusEvent : " + caseRoot_ + " " + caseName_;
    }

    public String caseRoot()
    {
        return caseRoot_;
    }

    public String caseName()
    {
        return caseName_;
    }

    public ICaseServer caseServer()
    {
        return caseServer_;
    }

    public ICasePostServer casePostServer()
    {
        return casePostServer_;
    }

    public ICaseBrowser caseBrowser()
    {
        return caseBrowser_;
    }
}



