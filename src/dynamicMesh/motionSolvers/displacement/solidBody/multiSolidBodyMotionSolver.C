/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "multiSolidBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "cellZoneMesh.H"
#include "boolList.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiSolidBodyMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        multiSolidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionSolver::multiSolidBodyMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    points0MotionSolver(mesh, dict, typeName)
{
    zoneIDs_.setSize(coeffDict().size());
    SBMFs_.setSize(coeffDict().size());
    pointIDs_.setSize(coeffDict().size());
    label zonei = 0;

    for (const entry& dEntry : coeffDict())
    {
        if (dEntry.isDict())
        {
            const word& zoneName = dEntry.keyword();
            const dictionary& subDict = dEntry.dict();

            zoneIDs_[zonei] = mesh.cellZones().findZoneID(zoneName);

            if (zoneIDs_[zonei] == -1)
            {
                FatalIOErrorInFunction(coeffDict())
                    << "Cannot find cellZone named " << zoneName
                    << ". Valid zones are "
                    << flatOutput(mesh.cellZones().names())
                    << exit(FatalIOError);
            }

            SBMFs_.set
            (
                zonei,
                solidBodyMotionFunction::New(subDict, mesh.time())
            );

            // Collect points of cell zone.
            const cellZone& cz = mesh.cellZones()[zoneIDs_[zonei]];

            boolList movePts(mesh.nPoints(), false);

            forAll(cz, i)
            {
                label celli = cz[i];
                const cell& c = mesh.cells()[celli];
                forAll(c, j)
                {
                    const face& f = mesh.faces()[c[j]];
                    forAll(f, k)
                    {
                        label pointi = f[k];
                        movePts[pointi] = true;
                    }
                }
            }

            syncTools::syncPointList(mesh, movePts, orEqOp<bool>(), false);

            DynamicList<label> ptIDs(mesh.nPoints());
            forAll(movePts, i)
            {
                if (movePts[i])
                {
                    ptIDs.append(i);
                }
            }

            pointIDs_[zonei].transfer(ptIDs);

            Info<< "Applying solid body motion " << SBMFs_[zonei].type()
                << " to "
                << returnReduce(pointIDs_[zonei].size(), sumOp<label>())
                << " points of cellZone " << zoneName << endl;

            zonei++;
        }
    }
    zoneIDs_.setSize(zonei);
    SBMFs_.setSize(zonei);
    pointIDs_.setSize(zonei);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiSolidBodyMotionSolver::~multiSolidBodyMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::multiSolidBodyMotionSolver::curPoints() const
{
    tmp<pointField> ttransformedPts(new pointField(mesh().points()));
    pointField& transformedPts = ttransformedPts.ref();

    forAll(zoneIDs_, i)
    {
        const labelList& zonePoints = pointIDs_[i];

        UIndirectList<point>(transformedPts, zonePoints) = transformPoints
        (
            SBMFs_[i].transformation(),
            pointField(points0_, zonePoints)
        );
    }

    return ttransformedPts;
}


// ************************************************************************* //
