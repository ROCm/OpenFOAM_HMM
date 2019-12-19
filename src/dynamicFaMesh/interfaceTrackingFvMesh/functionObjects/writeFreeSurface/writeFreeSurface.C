/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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

\*----------------------------------------------------------------------------*/

#include "writeFreeSurface.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"
#include "interfaceTrackingFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(writeFreeSurface, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeFreeSurface,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::writeFreeSurface::writeData()
{
    if (time_.outputTime())
    {
        const fvMesh& mesh =
            time_.lookupObject<fvMesh>(polyMesh::defaultRegion);

        interfaceTrackingFvMesh& itm =
            refCast<interfaceTrackingFvMesh>
            (
                const_cast<dynamicFvMesh&>
                (
                    mesh.lookupObject<dynamicFvMesh>("fvSolution")
                )
            );

        itm.writeVTKControlPoints();
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::writeFreeSurface::writeFreeSurface
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(runTime),
    regionName_(polyMesh::defaultRegion)
{
    Info<< "Creating " << this->name() << " function object." << endl;

    dict.readIfPresent("region", regionName_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::writeFreeSurface::start()
{
    return writeData();
}


bool Foam::writeFreeSurface::execute()
{
    return writeData();
}


bool Foam::writeFreeSurface::read(const dictionary& dict)
{
    dict.readIfPresent("region", regionName_);

    return true;
}


// ************************************************************************* //
