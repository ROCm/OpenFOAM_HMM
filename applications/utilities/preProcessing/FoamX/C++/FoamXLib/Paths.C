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

// Foam header files.
#include "OSspecific.H"

// FoamX header files.
#include "Paths.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Paths
{
    Foam::fileName config(Foam::getEnv("FOAMX_CONFIG"));
    //Foam::fileName tmp(Foam::dotFoam("tmp"));
    Foam::fileName tmp(config);

    Foam::fileName solversPath()
    {
        return Foam::getEnv("WM_PROJECT_DIR")/"applications/solvers";
    }
    Foam::fileName solvers(solversPath());

    Foam::fileName userSolversPath()
    {
        return Foam::getEnv("WM_PROJECT_USER_DIR")/"applications/solvers";
    }
    Foam::fileName userSolvers(userSolversPath());

    Foam::fileName utilitiesPath()
    {
        return Foam::getEnv("WM_PROJECT_DIR")/"applications/utilities";
    }
    Foam::fileName utilities(utilitiesPath());

    Foam::fileName userUtilitiesPath()
    {
        return Foam::getEnv("WM_PROJECT_USER_DIR")/"applications/utilities";
    }
    Foam::fileName userUtilities(userUtilitiesPath());
}


// ************************************************************************* //
