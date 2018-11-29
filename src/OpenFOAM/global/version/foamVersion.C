/*-------------------------------*- C++ -*-----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "foamVersion.H"
#include "messageStream.H"

// Static data members are constructed in global.Cver

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::foamVersion::patched()
{
    // Patch-level, when defined (non-zero) and not some @TOKEN@ rubbish
    return
    (
        foamVersion::patch.size() && foamVersion::patch[0] != '@'
     && (foamVersion::patch.size() > 1 || foamVersion::patch[0] != '0')
    );
}


void Foam::foamVersion::printBuildInfo(const bool full)
{
    // Can use #if OPENFOAM directly

    Info<< "Using: OpenFOAM-" << foamVersion::version.c_str()
        << " (" << OPENFOAM << ") (see www.OpenFOAM.com)" << nl
        << "Build: " << foamVersion::build.c_str();

    if (foamVersion::patched())
    {
        // Patch-level, when defined
        Info<< " (patch=" << foamVersion::patch.c_str() << ')';
    }
    Info<< nl;

    if (full)
    {
        Info<< "Arch:  " << foamVersion::buildArch.c_str() << nl;
    }
}


// ************************************************************************* //
