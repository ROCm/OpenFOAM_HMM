/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "dynamicMultiMotionSolverFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boolList.H"
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMultiMotionSolverFvMesh::dynamicMultiMotionSolverFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io)
{
    IOdictionary dynDict
    (
        IOobject
        (
            "dynamicMeshDict",
            io.time().constant(),
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
    label zoneI = 0;

    forAllConstIter(dictionary, dynamicMeshCoeffs, iter)
    {
        if (iter().isDict())
        {
            const dictionary& subDict = iter().dict();

            word zoneName(subDict.lookup("cellZone"));

            zoneIDs_[zoneI] = cellZones().findZoneID(zoneName);

            if (zoneIDs_[zoneI] == -1)
            {
                FatalIOErrorInFunction
                (
                    dynamicMeshCoeffs
                )   << "Cannot find cellZone named " << zoneName
                    << ". Valid zones are " << cellZones().names()
                    << exit(FatalIOError);
            }

            IOobject io(dynDict);
            io.readOpt() = IOobject::NO_READ;

            motionPtr_.set
            (
                zoneI,
                motionSolver::New
                (
                    *this,
                    IOdictionary(io, subDict)
                )
            );

            // Collect points of cell zone.
            const cellZone& cz = cellZones()[zoneIDs_[zoneI]];

            boolList movePts(nPoints(), false);

            forAll(cz, i)
            {
                label cellI = cz[i];
                const cell& c = cells()[cellI];
                forAll(c, j)
                {
                    const face& f = faces()[c[j]];
                    forAll(f, k)
                    {
                        label pointI = f[k];
                        movePts[pointI] = true;
                    }
                }
            }

            syncTools::syncPointList(*this, movePts, orEqOp<bool>(), false);

            DynamicList<label> ptIDs(nPoints());
            forAll(movePts, i)
            {
                if (movePts[i])
                {
                    ptIDs.append(i);
                }
            }

            pointIDs_[zoneI].transfer(ptIDs);

            Info<< "Applying motionSolver " << motionPtr_[zoneI].type()
                << " to "
                << returnReduce(pointIDs_[zoneI].size(), sumOp<label>())
                << " points of cellZone " << zoneName << endl;

            zoneI++;
        }
    }
    zoneIDs_.setSize(zoneI);
    motionPtr_.setSize(zoneI);
    pointIDs_.setSize(zoneI);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicMultiMotionSolverFvMesh::~dynamicMultiMotionSolverFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicMultiMotionSolverFvMesh::update()
{
    pointField transformedPts(points());

    forAll(motionPtr_, zoneI)
    {
        tmp<pointField> tnewPoints(motionPtr_[zoneI].newPoints());
        const pointField& newPoints = tnewPoints();

        const labelList& zonePoints = pointIDs_[zoneI];
        forAll(zonePoints, i)
        {
            label pointI = zonePoints[i];
            transformedPts[pointI] = newPoints[pointI];
        }
    }

    fvMesh::movePoints(transformedPts);

    static bool hasWarned = false;

    if (foundObject<volVectorField>("U"))
    {
        const_cast<volVectorField&>(lookupObject<volVectorField>("U"))
            .correctBoundaryConditions();
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
