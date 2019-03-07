/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "cyclicACMIPolyPatch.H"
#include "SubField.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicACMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicACMIPolyPatch, dictionary);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cyclicACMIPolyPatch::reportCoverage
(
    const word& name,
    const scalarField& weightSum
) const
{
    label nUncovered = 0;
    label nCovered = 0;
    for (const scalar sum : weightSum)
    {
        if (sum < tolerance_)
        {
            ++nUncovered;
        }
        else if (sum > scalar(1) - tolerance_)
        {
            ++nCovered;
        }
    }
    reduce(nUncovered, sumOp<label>());
    reduce(nCovered, sumOp<label>());
    label nTotal = returnReduce(weightSum.size(), sumOp<label>());

    Info<< "ACMI: Patch " << name << " uncovered/blended/covered = "
        << nUncovered << ", " << nTotal-nUncovered-nCovered
        << ", " << nCovered << endl;
}


void Foam::cyclicACMIPolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod
) const
{
    resetAMI(boundaryMesh().mesh().points(), AMIMethod);
}


void Foam::cyclicACMIPolyPatch::resetAMI
(
    const UList<point>& points,
    const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod
) const
{
    if (!owner())
    {
        return;
    }

    const polyPatch& nonOverlapPatch = this->nonOverlapPatch();

    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::resetAMI : recalculating weights"
            << " for " << name() << " and " << nonOverlapPatch.name()
            << endl;
    }


WarningInFunction<< "DEACTIVATED clearGeom()" << endl;
//    if (boundaryMesh().mesh().hasCellCentres())
    if (0 && boundaryMesh().mesh().hasCellCentres())
    {
        if (debug)
        {
            Pout<< "cyclicACMIPolyPatch::resetAMI : clearing cellCentres"
                << " for " << name() << " and " << nonOverlapPatch.name()
                << endl;
        }

        //WarningInFunction
        //    << "The mesh already has cellCentres calculated when"
        //    << " resetting ACMI " << name() << "." << endl
        //    << "This is a problem since ACMI adapts the face areas"
        //    << " (to close cells) so this has" << endl
        //    << "to be done before cell centre calculation." << endl
        //    << "This can happen if e.g. the cyclicACMI is after"
        //    << " any processor patches in the boundary." << endl;
        const_cast<polyMesh&>
        (
            boundaryMesh().mesh()
        ).primitiveMesh::clearGeom();
    }


    // Trigger re-building of faceAreas
    (void)boundaryMesh().mesh().faceAreas();


    // Calculate the AMI using partial face-area-weighted. This leaves
    // the weights as fractions of local areas (sum(weights) = 1 means
    // face is fully covered)
    cyclicAMIPolyPatch::resetAMI
    (
        points,
        AMIPatchToPatchInterpolation::imPartialFaceAreaWeight
    );

    const AMIPatchToPatchInterpolation& AMI = this->AMI();

    // Output some statistics
    reportCoverage("source", AMI.srcWeightsSum());
    reportCoverage("target", AMI.tgtWeightsSum());

    // Set the mask fields
    // Note:
    // - assumes that the non-overlap patches are decomposed using the same
    //   decomposition as the coupled patches (per side)
    srcMask_ = min(scalar(1), max(scalar(0), AMI.srcWeightsSum()));
    tgtMask_ = min(scalar(1), max(scalar(0), AMI.tgtWeightsSum()));

    if (debug)
    {
        Pout<< "resetAMI" << endl;
        {
            const cyclicACMIPolyPatch& patch = *this;
            Pout<< "patch:" << patch.name() << " size:" << patch.size()
                << " non-overlap patch: " << patch.nonOverlapPatch().name()
                << " size:" << patch.nonOverlapPatch().size()
                << " mask size:" << patch.srcMask().size() << endl;
        }
        {
            const cyclicACMIPolyPatch& patch = this->neighbPatch();
            Pout<< "patch:" << patch.name() << " size:" << patch.size()
                << " non-overlap patch: " << patch.nonOverlapPatch().name()
                << " size:" << patch.nonOverlapPatch().size()
                << " mask size:" << patch.neighbPatch().tgtMask().size()
                << endl;
        }
    }
}


void Foam::cyclicACMIPolyPatch::scalePatchFaceAreas()
{
    if (!owner() || !canResetAMI())
    {
        return;
    }

    scalePatchFaceAreas(*this);
    scalePatchFaceAreas(this->neighbPatch());
}


void Foam::cyclicACMIPolyPatch::scalePatchFaceAreas
(
    const cyclicACMIPolyPatch& acmipp
)
{
    // Primitive patch face areas have been cleared/reset based on the raw
    // points - need to reset to avoid double-accounting of face areas

    DebugPout
        << "rescaling non-overlap patch areas" << endl;

    const scalar maxTol = scalar(1) - tolerance_;
    const scalarField& mask = acmipp.mask();

    const polyPatch& nonOverlapPatch = acmipp.nonOverlapPatch();
    vectorField::subField noSf = nonOverlapPatch.faceAreas();

    if (mask.size() != noSf.size())
    {
        WarningInFunction
            << "Inconsistent sizes for patch: " << acmipp.name()
            << " - not manipulating patches" << nl
            << " - size: " << size() << nl
            << " - non-overlap patch size: " << noSf.size() << nl
            << " - mask size: " << mask.size() << nl
            << "This is OK for decomposition but should be considered fatal "
            << "at run-time" << endl;

        return;
    }

    forAll(noSf, facei)
    {
        const scalar w = min(maxTol, max(tolerance_, mask[facei]));
        noSf[facei] *= scalar(1) - w;
    }

    if (!createAMIFaces_)
    {
        // Note: for topological update (createAMIFaces_ = true)
        // AMI coupled patch face areas are updated as part of the topological
        // updates, e.g. by the calls to cyclicAMIPolyPatch's setTopology and
        // initMovePoints
        DebugPout
            << "scaling coupled patch areas" << endl;

        // Scale the coupled patch face areas
        vectorField::subField Sf = acmipp.faceAreas();

        forAll(Sf, facei)
        {
            Sf[facei] *= max(tolerance_, mask[facei]);
        }

        // Re-normalise the weights since the effect of overlap is already
        // accounted for in the area
        auto& weights = const_cast<scalarListList&>(acmipp.weights());
        auto& weightsSum = const_cast<scalarField&>(acmipp.weightsSum());
        forAll(weights, i)
        {
            scalarList& wghts = weights[i];
            if (wghts.size())
            {
                scalar& sum = weightsSum[i];

                forAll(wghts, j)
                {
                    wghts[j] /= sum;
                }
                sum = 1.0;
            }
        }
    }
}


void Foam::cyclicACMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::initGeometry : " << name() << endl;
    }

    // Note: calculates transformation and triggers face centre calculation
    cyclicAMIPolyPatch::initGeometry(pBufs);

    // On start-up there are no topological updates so scale the face areas
    // - Note: resetAMI called by cyclicAMIPolyPatch::initGeometry
    scalePatchFaceAreas();
}


void Foam::cyclicACMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::calcGeometry : " << name() << endl;
    }
    cyclicAMIPolyPatch::calcGeometry(pBufs);
}


void Foam::cyclicACMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::initMovePoints : " << name() << endl;
    }

    // Note: calculates transformation and triggers face centre calculation
    // - Note: resetAMI called by cyclicAMIPolyPatch::initMovePoints
    cyclicAMIPolyPatch::initMovePoints(pBufs, p);

    scalePatchFaceAreas();
}


void Foam::cyclicACMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::movePoints : " << name() << endl;
    }

    // When topology is changing, this will scale the duplicate AMI faces
    cyclicAMIPolyPatch::movePoints(pBufs, p);
}


void Foam::cyclicACMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::initUpdateMesh : " << name() << endl;
    }
    cyclicAMIPolyPatch::initUpdateMesh(pBufs);
}


void Foam::cyclicACMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::updateMesh : " << name() << endl;
    }
    cyclicAMIPolyPatch::updateMesh(pBufs);
}


void Foam::cyclicACMIPolyPatch::clearGeom()
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::clearGeom : " << name() << endl;
    }
    cyclicAMIPolyPatch::clearGeom();
}


const Foam::scalarField& Foam::cyclicACMIPolyPatch::srcMask() const
{
    return srcMask_;
}


const Foam::scalarField& Foam::cyclicACMIPolyPatch::tgtMask() const
{
    return tgtMask_;
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
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
    nonOverlapPatchName_(word::null),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_()
{
    AMIRequireMatch_ = false;

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch(name, dict, index, bm, patchType),
    nonOverlapPatchName_(dict.lookup("nonOverlapPatch")),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_()
{
    AMIRequireMatch_ = false;

    if (nonOverlapPatchName_ == name)
    {
        FatalIOErrorInFunction(dict)
            << "Non-overlapping patch name " << nonOverlapPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_()
{
    AMIRequireMatch_ = false;

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName,
    const word& nonOverlapPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    nonOverlapPatchName_(nonOverlapPatchName),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_()
{
    AMIRequireMatch_ = false;

    if (nonOverlapPatchName_ == name())
    {
        FatalErrorInFunction
            << "Non-overlapping patch name " << nonOverlapPatchName_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_()
{
    AMIRequireMatch_ = false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::cyclicACMIPolyPatch& Foam::cyclicACMIPolyPatch::neighbPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];
    return refCast<const cyclicACMIPolyPatch>(pp);
}


Foam::label Foam::cyclicACMIPolyPatch::nonOverlapPatchID() const
{
    if (nonOverlapPatchID_ == -1)
    {
        nonOverlapPatchID_ =
            this->boundaryMesh().findPatchID(nonOverlapPatchName_);

        if (nonOverlapPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal non-overlapping patch name " << nonOverlapPatchName_
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        if (nonOverlapPatchID_ < index())
        {
            FatalErrorInFunction
                << "Boundary ordering error: " << type()
                << " patch must be defined prior to its non-overlapping patch"
                << nl
                << type() << " patch: " << name() << ", ID:" << index() << nl
                << "Non-overlap patch: " << nonOverlapPatchName_
                << ", ID:" << nonOverlapPatchID_ << nl
                << exit(FatalError);
        }

        const polyPatch& noPp = this->boundaryMesh()[nonOverlapPatchID_];

        bool ok = true;

        if (size() == noPp.size())
        {
            const scalarField magSf(mag(faceAreas()));
            const scalarField noMagSf(mag(noPp.faceAreas()));

            forAll(magSf, facei)
            {
                scalar ratio = mag(magSf[facei]/(noMagSf[facei] + ROOTVSMALL));

                if (ratio - 1 > tolerance_)
                {
                    ok = false;
                    break;
                }
            }
        }
        else
        {
            ok = false;
        }

        if (!ok)
        {
            FatalErrorInFunction
                << "Inconsistent ACMI patches " << name() << " and "
                << noPp.name() << ".  Patches should have identical topology"
                << exit(FatalError);
        }
    }

    return nonOverlapPatchID_;
}


void Foam::cyclicACMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{
    cyclicAMIPolyPatch::initOrder(pBufs, pp);
}


bool Foam::cyclicACMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    return cyclicAMIPolyPatch::order(pBufs, pp, faceMap, rotation);
}


void Foam::cyclicACMIPolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);

    os.writeEntry("nonOverlapPatch", nonOverlapPatchName_);
}


// ************************************************************************* //
