/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2008-2011 OpenCFD Ltd.
    \\/      M anipulation   |
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

Description
    Initialises fields for a molecular dynamics (MD) simulation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "polyatomicCloud.H"
#include "monoatomicCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    IOdictionary mdInitialiseDict
    (
        IOobject
        (
            "mdInitialiseDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    word polyCloudName("polyatomicCloud");

    const dictionary& polyDict(mdInitialiseDict.subDict(polyCloudName));

    IOdictionary polyIdListDict
    (
        IOobject
        (
            polyCloudName + "_idList",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    potential polyPot
    (
        mesh,
        polyDict,
        IOdictionary
        (
            IOobject
            (
                polyCloudName +  "Properties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ),
        polyIdListDict
    );

    polyatomicCloud poly
    (
        polyCloudName,
        mesh,
        polyPot,
        polyDict
    );

    Info<< nl << returnReduce(poly.size(), sumOp<label>()) << " added to "
        << poly.name()
        << nl << endl;

    word monoCloudName("monoatomicCloud");

    const dictionary& monoDict(mdInitialiseDict.subDict(monoCloudName));

    IOdictionary monoIdListDict
    (
        IOobject
        (
            monoCloudName + "_idList",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    potential monoPot
    (
        mesh,
        monoDict,
        IOdictionary
        (
            IOobject
            (
                monoCloudName +  "Properties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ),
        monoIdListDict
    );

    monoatomicCloud mono
    (
        monoCloudName,
        mesh,
        monoPot,
        monoDict
    );

    Info<< nl << returnReduce(mono.size(), sumOp<label>()) << " added to "
        << mono.name()
        << nl << endl;

    if (!mesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing."
            << nl << exit(FatalError);
    }

    Info<< nl << "ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
