/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "dynamicMultiMotionSolverFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "bitSet.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMultiMotionSolverFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMultiMotionSolverFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMultiMotionSolverFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMultiMotionSolverFvMesh::dynamicMultiMotionSolverFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicFvMesh(io, doInit)
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}


bool Foam::dynamicMultiMotionSolverFvMesh::init(const bool doInit)
{
    if (doInit)
    {
        dynamicFvMesh::init(doInit);
    }

    IOdictionary dynDict
    (
        IOobject
        (
            "dynamicMeshDict",
            time().constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );
    const dictionary& dynamicMeshCoeffs = dynDict.subDict(typeName + "Coeffs");

    zoneIDs_.setSize(dynamicMeshCoeffs.size());
    motionPtr_.setSize(dynamicMeshCoeffs.size());
    pointIDs_.setSize(dynamicMeshCoeffs.size());

    label zonei = 0;

    bitSet movePts;

    for (const entry& dEntry : dynamicMeshCoeffs)
    {
        if (dEntry.isDict())
        {
            const dictionary& subDict = dEntry.dict();

            const word zoneName(subDict.get<word>("cellZone"));

            zoneIDs_[zonei] = cellZones().findZoneID(zoneName);

            if (zoneIDs_[zonei] == -1)
            {
                FatalIOErrorInFunction(dynamicMeshCoeffs)
                    << "Cannot find cellZone named " << zoneName
                    << ". Valid zones are " << cellZones().names()
                    << exit(FatalIOError);
            }

            IOobject io(dynDict);
            io.readOpt(IOobject::NO_READ);

            motionPtr_.set
            (
                zonei,
                motionSolver::New
                (
                    *this,
                    IOdictionary(io, subDict)
                )
            );


            // Collect points of cell zone.

            movePts.reset();
            movePts.resize(nPoints());

            for (const label celli : cellZones()[zoneIDs_[zonei]])
            {
                for (const label facei : cells()[celli])
                {
                    movePts.set(faces()[facei]);  // set multiple points
                }
            }

            syncTools::syncPointList
            (
                *this, movePts, orEqOp<unsigned int>(), 0u
            );

            pointIDs_[zonei] = movePts.sortedToc();

            Info<< "Applying motionSolver " << motionPtr_[zonei].type()
                << " to "
                << returnReduce(pointIDs_[zonei].size(), sumOp<label>())
                << " points of cellZone " << zoneName << endl;

            ++zonei;
        }
    }
    zoneIDs_.setSize(zonei);
    motionPtr_.setSize(zonei);
    pointIDs_.setSize(zonei);

    // Assume changed ...
    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicMultiMotionSolverFvMesh::update()
{
    pointField transformedPts(points());

    forAll(motionPtr_, zonei)
    {
        tmp<pointField> tnewPoints(motionPtr_[zonei].newPoints());
        const pointField& newPoints = tnewPoints();

        for (const label pointi : pointIDs_[zonei])
        {
            transformedPts[pointi] = newPoints[pointi];
        }
    }

    fvMesh::movePoints(transformedPts);

    static bool hasWarned = false;

    volVectorField* Uptr = getObjectPtr<volVectorField>("U");

    if (Uptr)
    {
        Uptr->correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningInFunction
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }

    return true;
}


// ************************************************************************* //
