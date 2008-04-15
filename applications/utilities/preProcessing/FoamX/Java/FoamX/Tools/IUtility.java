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
package FoamX.Tools;

import FoamXServer.CaseBrowser.ICaseBrowser;
import FoamXServer.CaseServer.ICaseServer;
import FoamXServer.UtilityDescriptor;

public interface IUtility
{
    public void prepareUtility();

    // Utilities called at browser level are invoked with
    // caseBrowser instance and -if possible- with root and casename
    public void invokeBrowserUtility
    (
        UtilityDescriptor toolDesc,
        ICaseBrowser caseBrowser,
        String caseRoot,
        String caseName
    );

    // Utilities called at case editor level are invoked with
    // caseServer instance
    public void invokeCaseUtility
    (
        UtilityDescriptor toolDesc,
        ICaseServer caseServer,
        String caseRoot,
        String caseName
    );
}

