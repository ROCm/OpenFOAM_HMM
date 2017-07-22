/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

Application
    foamHelp

Group
    grpMiscUtilities

Description
    Top level wrapper utility around foam help utilities

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "helpType.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "addToolOption.H"

    // Intercept request for help
    if ((argc > 1) && (strcmp(argv[1], "-help") == 0))
    {
        #include "setRootCase.H"
    }

    if (argc < 2)
    {
        FatalError
            << "No help utility has been supplied" << nl
            << exit(FatalError);
    }

    word utilityName = argv[1];
    Foam::autoPtr<Foam::helpType> utility
    (
        helpType::New(utilityName)
    );

    utility().init();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    utility().execute(args, mesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
