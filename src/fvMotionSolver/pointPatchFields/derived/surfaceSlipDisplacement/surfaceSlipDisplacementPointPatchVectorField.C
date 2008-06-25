/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "surfaceSlipDisplacementPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "transformField.H"
#include "fvMesh.H"
#include "displacementLaplacianFvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char*
NamedEnum<surfaceSlipDisplacementPointPatchVectorField::followMode, 3>::
names[] =
{
    "nearest",
    "pointNormal",
    "fixedNormal"
};

const NamedEnum<surfaceSlipDisplacementPointPatchVectorField::followMode, 3>
    surfaceSlipDisplacementPointPatchVectorField::followModeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    pointPatchVectorField(p, iF),
    surfaceNames_(),
    projectMode_(NEAREST),
    projectDir_(vector::zero),
    wedgePlane_(-1)
{}


surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    pointPatchVectorField(p, iF, dict),
    surfaceNames_(dict.lookup("projectSurfaces")),
    projectMode_(followModeNames_.read(dict.lookup("followMode"))),
    projectDir_(dict.lookup("projectDirection")),
    wedgePlane_(readLabel(dict.lookup("wedgePlane"))),
    frozenPointsZone_(dict.lookup("frozenPointsZone"))
{}


surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const surfaceSlipDisplacementPointPatchVectorField& ppf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper&
)
:
    pointPatchVectorField(p, iF),
    surfaceNames_(ppf.surfaceNames()),
    projectMode_(ppf.projectMode()),
    projectDir_(ppf.projectDir()),
    wedgePlane_(ppf.wedgePlane()),
    frozenPointsZone_(ppf.frozenPointsZone())
{}


surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const surfaceSlipDisplacementPointPatchVectorField& ppf
)
:
    pointPatchVectorField(ppf),
    surfaceNames_(ppf.surfaceNames()),
    projectMode_(ppf.projectMode()),
    projectDir_(ppf.projectDir()),
    wedgePlane_(ppf.wedgePlane()),
    frozenPointsZone_(ppf.frozenPointsZone())
{}


surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const surfaceSlipDisplacementPointPatchVectorField& ppf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    pointPatchVectorField(ppf, iF),
    surfaceNames_(ppf.surfaceNames()),
    projectMode_(ppf.projectMode()),
    projectDir_(ppf.projectDir()),
    wedgePlane_(ppf.wedgePlane()),
    frozenPointsZone_(ppf.frozenPointsZone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const triSurfaceMeshes& surfaceSlipDisplacementPointPatchVectorField::
surfaces() const
{
    if (!surfacesPtr_.valid())
    {
        surfacesPtr_.reset
        (
            new triSurfaceMeshes
            (
                IOobject
                (
                    "abc",                              // dummy name
                    db().time().constant(),             // directory
                    "triSurface",                       // instance
                    db().time(),                        // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfaceNames_
            )
        );
    }
    return surfacesPtr_();
}


void surfaceSlipDisplacementPointPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    const polyMesh& mesh = patch().boundaryMesh().mesh()();

    //const scalar deltaT = mesh.time().deltaT().value();

    // Construct large enough vector in direction of projectDir so
    // we're guaranteed to hit something.

    const scalar projectLen = mag(mesh.bounds().max()-mesh.bounds().min());

    // For case of fixed projection vector:
    vector projectVec;
    if (projectMode_ == FIXEDNORMAL)
    {
        vector n = projectDir_/mag(projectDir_);
        projectVec = projectLen*n;
    }

    //- Per point projection vector:

    const pointField& localPoints = patch().localPoints();
    const labelList& meshPoints = patch().meshPoints();

    vectorField displacement(this->patchInternalField());


    // Get fixed points (bit of a hack)
    const pointZone* zonePtr = NULL;

    if (frozenPointsZone_.size() > 0)
    {
        const pointZoneMesh& pZones = mesh.pointZones();

        zonePtr = &pZones[pZones.findZoneID(frozenPointsZone_)];

        Pout<< "surfaceSlipDisplacementPointPatchVectorField : Fixing all "
            << zonePtr->size() << " points in pointZone " << zonePtr->name()
            << endl;
    }

    // Get the starting locations from the motionSolver
    const displacementLaplacianFvMotionSolver& motionSolver =
        mesh.lookupObject<displacementLaplacianFvMotionSolver>
        (
            "dynamicMeshDict"
        );
    const pointField& points0 = motionSolver.points0();

    forAll(localPoints, i)
    {
        if (zonePtr && (zonePtr->whichPoint(meshPoints[i]) >= 0))
        {
            // Fixed point. Reset to point0 location.

            //Pout<< "    Fixed point:" << meshPoints[i]
            //    << " coord:" << localPoints[i]
            //    << " should be at:" << points0[meshPoints[i]]
            //    << endl;
            displacement[i] = points0[meshPoints[i]] - localPoints[i];
        }
        else
        {
            point start(points0[meshPoints[i]] + displacement[i]);

            scalar offset = 0;
            pointIndexHit intersection;

            if (projectMode_ == NEAREST)
            {
                surfaces().findNearest(start, sqr(projectLen), intersection);
            }
            else
            {
                // Check if already on surface
                surfaces().findNearest(start, sqr(SMALL), intersection);

                if (!intersection.hit())
                {
                    // No nearest found. Do intersection

                    if (projectMode_ == POINTNORMAL)
                    {
                        projectVec = projectLen*patch().pointNormals()[i];
                    }

                    // Knock out any wedge component
                    if (wedgePlane_ >= 0 && wedgePlane_ <= vector::nComponents)
                    {
                        offset = start[wedgePlane_];
                        start[wedgePlane_] = 0;
                        projectVec[wedgePlane_] = 0;
                    }

                    label rightSurf0, rightSurf1;
                    pointIndexHit rightHit0, rightHit1;
                    surfaces().findNearestIntersection
                    (
                        start,
                        start+projectVec,
                        rightSurf0,
                        rightHit0,
                        rightSurf1,
                        rightHit1
                    );

                    // Do intersection
                    label leftSurf0, leftSurf1;
                    pointIndexHit leftHit0, leftHit1;
                    surfaces().findNearestIntersection
                    (
                        start,
                        start-projectVec,
                        leftSurf0,
                        leftHit0,
                        leftSurf1,
                        leftHit1
                    );

                    if (rightHit0.hit())
                    {
                        if (leftHit0.hit())
                        {
                            if
                            (
                                magSqr(rightHit0.hitPoint()-start)
                              < magSqr(leftHit0.hitPoint()-start)
                            )
                            {
                                intersection = rightHit0;
                            }
                            else
                            {
                                intersection = leftHit0;
                            }
                        }
                        else
                        {
                            intersection = rightHit0;
                        }
                    }
                    else
                    {
                        if (leftHit0.hit())
                        {
                            intersection = leftHit0;
                        }
                    }
                }
            }

            // Update *this from intersection point

            if (intersection.hit())
            {
                point interPt = intersection.hitPoint();

                if (wedgePlane_ >= 0 && wedgePlane_ <= vector::nComponents)
                {
                    interPt[wedgePlane_] += offset;
                }
                displacement[i] = interPt-points0[meshPoints[i]];
            }
            else
            {
                Pout<< "    point:" << meshPoints[i]
                    << " coord:" << localPoints[i]
                    << "  did not find any intersection between ray from "
                    << start-projectVec << " to " << start+projectVec
                    << endl;
            }  
        }
    }

    // Get internal field to insert values into
    Field<vector>& iF = const_cast<Field<vector>&>(this->internalField());

    //setInInternalField(iF, motionU);
    setInInternalField(iF, displacement);

    pointPatchVectorField::evaluate(commsType);
}


void surfaceSlipDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchVectorField::write(os);
    os.writeKeyword("projectSurfaces") << surfaceNames_
        << token::END_STATEMENT << nl;
    os.writeKeyword("followMode") << followModeNames_[projectMode_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("projectDirection") << projectDir_
        << token::END_STATEMENT << nl;
    os.writeKeyword("wedgePlane") << wedgePlane_
        << token::END_STATEMENT << nl;
    os.writeKeyword("frozenPointsZone") << frozenPointsZone_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    surfaceSlipDisplacementPointPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
