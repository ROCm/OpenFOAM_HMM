/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "faMesh.H"
#include "polyMesh.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

bool Foam::faMesh::hasSystemFiles(const polyMesh& pMesh)
{
    // Expect
    // - system/faSchemes
    // - system/faSolution

    const fileOperation& fp = Foam::fileHandler();

    bool looksValid = true;

    // Global files: system/{faSchemes,faSolution}
    for
    (
        const word& expect
      : List<word>
        ({
            {"faSchemes"},
            {"faSolution"}
        })
    )
    {
        fileName found
        (
            fp.filePath
            (
                true,  // global
                IOobject
                (
                    expect,
                    pMesh.time().system(),
                    pMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                expect // typeName (ununsed?)
            )
        );

        if (found.empty())
        {
            looksValid = false;
        }
    }

    // Only needed on master
    Pstream::broadcast(looksValid);

    return looksValid;
}


bool Foam::faMesh::hasFiles(const polyMesh& pMesh)
{
    // As well as system/{faSchemes,faSolution}
    //
    // expect these:
    // - instance/faMesh/faceLabels
    // - instance/faMesh/faBoundary

    bool looksValid = hasSystemFiles(pMesh);

    if (looksValid)
    {
        const fileOperation& fp = Foam::fileHandler();

        // The geometry instance for faMesh/faceLabels
        // Must use READ_IF_PRESENT to avoid aborting if not available

        const word instance = pMesh.time().findInstance
        (
            pMesh.dbDir()/faMesh::meshSubDir,
            "faceLabels",
            IOobject::READ_IF_PRESENT
        );

        for
        (
            const wordPair& expect
          : List<wordPair>
            ({
                {"faceLabels", "labelList"},
                {"faBoundary", "faBoundaryMesh"}
            })
        )
        {
            const word& dataFile = expect.first();
            const word& dataClass = expect.second();

            fileName found
            (
                fp.filePath
                (
                    false,      // non-global
                    IOobject
                    (
                        dataFile,
                        instance,
                        faMesh::meshSubDir,
                        pMesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    dataClass   // typeName (ununsed?)
                )
            );

            if (found.empty())
            {
                looksValid = false;
            }
        }

        // Everybody needs it, or they all fail
        Pstream::reduceAnd(looksValid);
    }

    return looksValid;
}


Foam::autoPtr<Foam::faMesh> Foam::faMesh::TryNew(const polyMesh& pMesh)
{
    if (faMesh::hasFiles(pMesh))
    {
        return autoPtr<faMesh>::New(pMesh);
    }

    return nullptr;
}


// ************************************************************************* //
