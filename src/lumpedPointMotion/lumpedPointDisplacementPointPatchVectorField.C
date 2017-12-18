/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "lumpedPointDisplacementPointPatchVectorField.H"
#include "lumpedPointMovement.H"
#include "lumpedPointIOMovement.H"
#include "addToRunTimeSelectionTable.H"
#include "pointFields.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "Time.H"
#include "polyMesh.H"
#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        lumpedPointDisplacementPointPatchVectorField
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::labelList
Foam::lumpedPointDisplacementPointPatchVectorField::patchIds
(
    const pointVectorField& pvf
)
{
    const pointVectorField::Boundary& bf = pvf.boundaryField();

    DynamicList<label> patchLst(bf.size());
    forAll(bf, patchi)
    {
        // All patches of this type
        if (isA<patchType>(bf[patchi]))
        {
            patchLst.append(patchi);
            // or patchLst.append(bf[patchi].patch().index());
        }
    }

    return patchLst.shrink();
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

const Foam::pointField&
Foam::lumpedPointDisplacementPointPatchVectorField::points0() const
{
    const objectRegistry& obr = this->patch().boundaryMesh().mesh().db();

    // Obtain starting locations from the motionSolver
    return obr.lookupObject<displacementMotionSolver>
    (
        "dynamicMeshDict"
    ).points0();
}


const Foam::lumpedPointMovement&
Foam::lumpedPointDisplacementPointPatchVectorField::movement() const
{
    const objectRegistry& obr = this->patch().boundaryMesh().mesh().db();
    const lumpedPointIOMovement* ptr =
        lumpedPointIOMovement::lookupInRegistry(obr);

    if (ptr)
    {
        return *ptr; // Already exists
    }

    // create and register with this patch as the owner
    autoPtr<lumpedPointIOMovement> obj = lumpedPointIOMovement::New
    (
        obr,
        this->patch().index()
    );

    return objectRegistry::store(obj);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedPointDisplacementPointPatchVectorField::
lumpedPointDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF)
{}


Foam::lumpedPointDisplacementPointPatchVectorField::
lumpedPointDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict)
{}


Foam::lumpedPointDisplacementPointPatchVectorField::
lumpedPointDisplacementPointPatchVectorField
(
    const lumpedPointDisplacementPointPatchVectorField& pf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(pf, p, iF, mapper)
{}


Foam::lumpedPointDisplacementPointPatchVectorField::
lumpedPointDisplacementPointPatchVectorField
(
    const lumpedPointDisplacementPointPatchVectorField& pf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(pf, iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lumpedPointDisplacementPointPatchVectorField::
~lumpedPointDisplacementPointPatchVectorField()
{
    // de-register movement if in use and managed by this patch
    const lumpedPointIOMovement* ptr = lumpedPointIOMovement::lookupInRegistry
    (
        this->patch().boundaryMesh().mesh().db()
    );

    if (ptr && ptr->ownerId() == this->patch().index())
    {
        movement().coupler().shutdown();

        const_cast<lumpedPointIOMovement*>(ptr)->checkOut();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedPointDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    enum Time::stopAtControls action = Time::stopAtControls::saUnknown;

    const bool masterPatch = (movement().ownerId() == this->patch().index());
    if (masterPatch)
    {
        if (lumpedPointIOMovement::debug)
        {
            Pout<<"masterPatch: " << this->patch().index() << endl;
        }

        const polyMesh& mesh = this->patch().boundaryMesh().mesh().mesh();

        // need face 'zones' for calculating forces
        // likely need bounding box for the movement
        // -> do both now if required
        if (!movement().hasMapping())
        {
            const_cast<lumpedPointMovement&>(movement()).setMapping
            (
                mesh,
                // All patches of this type
                patchIds
                (
                    static_cast<const pointVectorField&>
                    (
                        this->internalField()
                    )
                ),
                this->points0()
            );
        }


        if
        (
            movement().coupler().initialized()
         || !movement().coupler().slaveFirst()
        )
        {
            // Synchronized for all processes
            List<vector> forces, moments;
            movement().forcesAndMoments(mesh, forces, moments);

            if (lumpedPointIOMovement::debug)
            {
                Pout<<"gatherForces: " << forces << " called from patch "
                    << this->patch().index() << endl;

                if (Pstream::master())
                {
                    Pout<<"output forces to file: "
                        << movement().locations() << " called from patch "
                        << this->patch().index() << nl
                        <<"# " << forces.size() << " force entries" << nl
                        <<"# fx fy fz" << nl
                        <<"output forces to file: "
                        << forces << " called from patch "
                        << this->patch().index() << endl;
                }
            }

            if (Pstream::master())
            {
                movement().writeData(forces, moments);

                // Signal external source to execute
                movement().coupler().useSlave();
            }
        }

        // Wait for slave to provide data (includes MPI barrier)
        // and catch any abort information sent from slave
        action = movement().coupler().waitForSlave();

        // Read data passed back from external source - includes MPI barrier
        const_cast<lumpedPointMovement&>(movement()).readState();
    }

    tmp<pointField> tdisp = movement().displacePoints
    (
        this->points0(),
        this->patch().meshPoints()
    );

    this->operator==(tdisp);

    fixedValuePointPatchField<vector>::updateCoeffs();

    // Process any abort information sent from slave
    if
    (
        action != this->db().time().stopAt()
     && action != Time::stopAtControls::saUnknown
    )
    {
        this->db().time().stopAt(action);
    }
}


void Foam::lumpedPointDisplacementPointPatchVectorField::write(Ostream& os)
const
{
    pointPatchField<vector>::write(os);
    writeEntry("value", os);
}


// ************************************************************************* //
