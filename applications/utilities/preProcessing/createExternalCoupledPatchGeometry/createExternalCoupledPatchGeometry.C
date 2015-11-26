/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

Application
    createExternalCoupledPatchGeometry.

Description
    Application to generate the patch geometry (points and faces) for use
    with the externalCoupled functionObject.

Usage
    - createExternalCoupledPatchGeometry \<patchGroup\> [OPTION]

    \param -commsDir \<commsDir\> \n
    Specify an alternative communications directory (default is comms
    in the case directory)

    \param -region \<name\> \n
    Specify an alternative mesh region.

    On execution, the combined patch geometry (points and faces) are output
    to the communications directory.

Note:
    The addressing is patch-local, i.e. point indices for each patch point
    used for face addressing starts at index 0.

SeeAlso
    externalCoupledFunctionObject

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "externalCoupledFunctionObject.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    argList::validArgs.append("patchGroup");
    argList::addOption
    (
        "commsDir",
        "dir",
        "specify alternate communications directory. default is 'comms'"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const wordRe patchGroup(args[1]);

    fileName commsDir(runTime.path()/"comms");
    args.optionReadIfPresent("commsDir", commsDir);

    externalCoupledFunctionObject::writeGeometry(mesh, commsDir, patchGroup);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
