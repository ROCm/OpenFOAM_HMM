/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "globalIndexAndTransform.H"
#include "coupledPolyPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * Private Static Data Members * * * * * * * * * * * //

const Foam::label Foam::globalIndexAndTransform::base_ = 32;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::globalIndexAndTransform::determineTransforms()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    transforms_ = List<vector>(6, vector::zero);

    label nextTrans = 0;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<coupledPolyPatch>(pp))
        {
            const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>(pp);

            if (cpp.separated())
            {
                const vectorField& sepVecs = cpp.separation();

                forAll(sepVecs, sVI)
                {
                    const vector& sepVec = sepVecs[sVI];

                    if (mag(sepVec) > SMALL)
                    {
                        scalar tol = coupledPolyPatch::matchTol*mag(sepVec);

                        if (min(mag(transforms_ - sepVec)) > tol)
                        {
                            transforms_[nextTrans++] = sepVec;
                        }

                        if (nextTrans > 6)
                        {
                            FatalErrorIn
                            (
                                 "void Foam::globalIndexAndTransform::"
                                 "determineTransforms()"
                            )
                                << "More than six unsigned transforms detected:"
                                << nl << transforms_
                                << exit(FatalError);
                        }
                    }
                }
            }
            else if (!cpp.parallel())
            {
                WarningIn
                (
                    "void Foam::globalIndexAndTransform::determineTransforms()"
                )
                    << "Non-parallel patches not supported." << endl;

                // Pout<< cpp.forwardT() << nl
                //     << cpp.reverseT() << nl
                //     << (cpp.forwardT() & cpp.reverseT())
                //     << endl;
            }
        }
    }

    List<List<vector> > allTransforms(Pstream::nProcs());

    allTransforms[Pstream::myProcNo()] = transforms_;

    Pstream::gatherList(allTransforms);

    if (Pstream::master())
    {
        transforms_ = List<vector>(3, vector::zero);

        label nextTrans = 0;

        forAll(allTransforms, procI)
        {
            const List<vector>& procTransVecs = allTransforms[procI];

            forAll(procTransVecs, pSVI)
            {
                const vector& procTransVec = procTransVecs[pSVI];

                if (mag(procTransVec) > SMALL)
                {
                    scalar tol = coupledPolyPatch::matchTol*mag(procTransVec);

                    if
                    (
                        min(mag(procTransVec - transforms_)) > tol
                     && min(mag(procTransVec + transforms_)) > tol
                    )
                    {
                        transforms_[nextTrans++] = procTransVec;
                    }

                    if (nextTrans > 3)
                    {
                        FatalErrorIn
                        (
                            "List<vector> allTransforms(const polyMesh& mesh)"
                        )
                            << "More than three independent basic "
                            << "transforms detected:" << nl
                            << allTransforms
                            << transforms_
                            << exit(FatalError);
                    }
                }
            }
        }

        transforms_.setSize(nextTrans);
    }

    Pstream::scatter(transforms_);
}


void Foam::globalIndexAndTransform::determineTransformPermutations()
{
    label nTransformPermutations = pow(3, transforms_.size());

    transformPermutations_.setSize(nTransformPermutations);

    forAll(transformPermutations_, tPI)
    {
        vector trans = vector::zero;

        label transformIndex = tPI;

        // Invert the ternary index encoding using repeated division by
        // three

        for (label b = 0; b < transforms_.size(); b++)
        {
            label w = (transformIndex % 3) - 1;

            transformIndex /= 3;

            trans += w*transforms_[b];
        }

        transformPermutations_[tPI] = trans;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalIndexAndTransform::globalIndexAndTransform
(
    const polyMesh& mesh
)
:
    mesh_(mesh)
{
    determineTransforms();

    determineTransformPermutations();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::globalIndexAndTransform::~globalIndexAndTransform()
{}


// ************************************************************************* //
