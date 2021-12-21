/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    writeActiveDesignVariables

Description
    Writes the active design variables based on the selected parameterisation
    scheme, as read from the meshMovement part of optimisationDict.
    Keeps a back-up of the original optimisationDict in
    system/optimisationDict.org, as comments in the dict will be lost.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "optMeshMovement.H"
#include "updateMethod.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary optDict
    (
        IOobject
        (
            "optimisationDict",
            mesh.time().caseSystem(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Back-up old optimisationDict, to maintain potential comments in it
    if (Pstream::master())
    {
        Foam::cp(optDict.objectPath(), optDict.objectPath() + ".org");
    }

    // Construct mesh movement object and grab active design variables
    // Will exit with error if not implemented for this type
    const dictionary& movementDict =
        optDict.subDict("optimisation").subDict("meshMovement");
    // Empty patch list will do
    labelList patchIDs(0);
    autoPtr<optMeshMovement> movementPtr
        (optMeshMovement::New(mesh, movementDict, patchIDs));
    const labelList activeDesignVariables =
        movementPtr().getActiveDesignVariables();

    // Construct update method to grab the type
    dictionary& updateMethodDict =
        optDict.subDict("optimisation").subDict("updateMethod");
    autoPtr<updateMethod> updMethod(updateMethod::New(mesh, updateMethodDict));

    // Add to appropriate dictionary (creates it if it does not exist)
    dictionary& coeffsDict = updateMethodDict.subDictOrAdd(updMethod().type());
    coeffsDict.add<labelList>
    (
        "activeDesignVariables",
        activeDesignVariables,
        true
    );

    // Write modified dictionary
    optDict.regIOobject::writeObject
    (
        IOstreamOption(IOstream::ASCII),
        true
    );

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
