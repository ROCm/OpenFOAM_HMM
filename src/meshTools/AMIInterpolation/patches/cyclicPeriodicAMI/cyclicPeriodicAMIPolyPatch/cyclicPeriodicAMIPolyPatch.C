/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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

#include "cyclicPeriodicAMIPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicPeriodicAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicPeriodicAMIPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        cyclicPeriodicAMIPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::cyclicPeriodicAMIPolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod
) const
{
    if (owner())
    {
        // Get the periodic patch
        const coupledPolyPatch& periodicPatch
        (
            refCast<const coupledPolyPatch>
            (
                boundaryMesh()[periodicPatchID()]
            )
        );

        // Create copies of both patches' points, transformed to the owner
        pointField thisPoints0(localPoints());
        pointField nbrPoints0(neighbPatch().localPoints());
        transformPosition(nbrPoints0);

        // Reset the stored number of periodic transformations to a lower
        // absolute value if possible
        if (nSectors_ > 0)
        {
            if (nTransforms_ > nSectors_/2)
            {
                nTransforms_ -= nSectors_;
            }
            else if (nTransforms_ < - nSectors_/2)
            {
                nTransforms_ += nSectors_;
            }
        }

        // Apply the stored number of periodic transforms
        for (label i = 0; i < nTransforms_; ++ i)
        {
            periodicPatch.transformPosition(thisPoints0);
        }
        for (label i = 0; i > nTransforms_; -- i)
        {
            periodicPatch.transformPosition(nbrPoints0);
        }

        // Create another copy
        pointField thisPoints(thisPoints0);
        pointField nbrPoints(nbrPoints0);

        // Create patches for all the points
        primitivePatch thisPatch0
        (
            SubList<face>(localFaces(), size()),
            thisPoints0
        );
        primitivePatch thisPatch
        (
            SubList<face>(localFaces(), size()),
            thisPoints
        );
        primitivePatch nbrPatch0
        (
            SubList<face>(neighbPatch().localFaces(), neighbPatch().size()),
            nbrPoints0
        );
        primitivePatch nbrPatch
        (
            SubList<face>(neighbPatch().localFaces(), neighbPatch().size()),
            nbrPoints
        );

        // Construct a new AMI interpolation between the initial patch locations
        AMIPtr_.reset
        (
            new AMIPatchToPatchInterpolation
            (
                thisPatch0,
                nbrPatch0,
                surfPtr(),
                faceAreaIntersect::tmMesh,
                false,
                AMIPatchToPatchInterpolation::imPartialFaceAreaWeight,
                AMILowWeightCorrection_,
                AMIReverse_
            )
        );

        // Number of geometry replications
        label iter(0);
        label nTransformsOld(nTransforms_);

        // Weight sum averages
        scalar srcSum(gAverage(AMIPtr_->srcWeightsSum()));
        scalar tgtSum(gAverage(AMIPtr_->tgtWeightsSum()));

        // Direction of geometry replication
        bool direction = nTransforms_ >= 0;

        // Increase in the source weight sum for the last iteration in the
        // opposite direction. If the current increase is less than this, the
        // direction is reversed.
        scalar srcSumDiff = 0;

        // Loop, replicating the geometry
        while
        (
            (iter < maxIter_)
         && (
                (1 - srcSum > matchTolerance())
             || (1 - tgtSum > matchTolerance())
            )
        )
        {
            if (direction)
            {
                periodicPatch.transformPosition(thisPoints);

                thisPatch.movePoints(thisPoints);

                AMIPtr_->append(thisPatch, nbrPatch0);
            }
            else
            {
                periodicPatch.transformPosition(nbrPoints);

                nbrPatch.movePoints(nbrPoints);

                AMIPtr_->append(thisPatch0, nbrPatch);
            }

            const scalar srcSumNew = gAverage(AMIPtr_->srcWeightsSum());
            const scalar srcSumDiffNew = srcSumNew - srcSum;

            if (srcSumDiffNew < srcSumDiff || srcSumDiffNew < SMALL)
            {
                direction = !direction;

                srcSumDiff = srcSumDiffNew;
            }

            srcSum = srcSumNew;
            tgtSum = gAverage(AMIPtr_->tgtWeightsSum());

            nTransforms_ += direction ? +1 : -1;

            ++ iter;
        }

        // Average the number of transformstions
        nTransforms_ = (nTransforms_ + nTransformsOld)/2;

        // Check that the match is complete
        if (iter == maxIter_)
        {
            FatalErrorIn
            (
                "void Foam::cyclicPeriodicAMIPolyPatch::resetPeriodicAMI"
                "("
                    "const AMIPatchToPatchInterpolation::interpolationMethod&"
                ") const"
            )
                << "Patches " << name() << " and " << neighbPatch().name()
                << " do not couple to within a tolerance of "
                << matchTolerance()
                << " when transformed according to the periodic patch "
                << periodicPatch.name() << "." << exit(FatalError);
        }

        // Check that at least one patch has a weight sum of one
        if
        (
            mag(1 - srcSum) > matchTolerance()
         && mag(1 - tgtSum) > matchTolerance()
        )
        {
            FatalErrorIn
            (
                "void Foam::cyclicPeriodicAMIPolyPatch::resetPeriodicAMI"
                "("
                    "const AMIPatchToPatchInterpolation::interpolationMethod&"
                ") const"
            )
                << "Patches " << name() << " and " << neighbPatch().name()
                << " do not overlap an integer number of times when transformed"
                << " according to the periodic patch "
                << periodicPatch.name() << "." << exit(FatalError);
        }

        // Normalise the weights
        AMIPtr_->normaliseWeights(true, false);
    }
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform
)
:
    cyclicAMIPolyPatch(name, size, start, index, bm, patchType, transform),
    periodicPatchName_(word::null),
    periodicPatchID_(-1),
    nTransforms_(0),
    nSectors_(0),
    maxIter_(36)
{}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch(name, dict, index, bm, patchType),
    periodicPatchName_(dict.lookup("periodicPatch")),
    periodicPatchID_(-1),
    nTransforms_(dict.lookupOrDefault<label>("nTransforms", 0)),
    nSectors_(dict.lookupOrDefault<label>("nSectors", 0)),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 36))
{}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const cyclicPeriodicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
    nTransforms_(pp.nTransforms_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_)
{}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const cyclicPeriodicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
    nTransforms_(pp.nTransforms_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_)
{}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const cyclicPeriodicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
    nTransforms_(pp.nTransforms_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicPeriodicAMIPolyPatch::~cyclicPeriodicAMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cyclicPeriodicAMIPolyPatch::periodicPatchID() const
{
    if (periodicPatchName_ == word::null)
    {
        periodicPatchID_ = -1;

        return periodicPatchID_;
    }

    if (periodicPatchID_ == -1)
    {
        periodicPatchID_ = this->boundaryMesh().findPatchID(periodicPatchName_);

        if (periodicPatchID_ == -1)
        {
            FatalErrorIn("cyclicPolyAMIPatch::periodicPatchID() const")
                << "Illegal periodicPatch name " << periodicPatchName_
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a coupled patch
        refCast<const coupledPolyPatch>
        (
            this->boundaryMesh()[periodicPatchID_]
        );
    }

    return periodicPatchID_;
}


void Foam::cyclicPeriodicAMIPolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);

    os.writeKeyword("periodicPatch") << periodicPatchName_
        << token::END_STATEMENT << nl;

    if (nTransforms_ != 0)
    {
        os.writeKeyword("nTransforms") << nTransforms_ <<
            token::END_STATEMENT << nl;
    }

    if (nSectors_ != 0)
    {
        os.writeKeyword("nSectors") << nSectors_ <<
            token::END_STATEMENT << nl;
    }

    if (maxIter_ != 36)
    {
        os.writeKeyword("maxIter") << maxIter_ << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
