/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "cyclicAMIPolyPatch.H"
#include "transformField.H"
#include "SubField.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "faceAreaIntersect.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, dictionary);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::label Foam::cyclicAMIPolyPatch::findFaceMaxRadius
(
    const pointField& faceCentres
) const
{
    // Determine a face furthest away from the axis

    const scalarField magRadSqr
    (
        magSqr((faceCentres - rotationCentre_) ^ rotationAxis_)
    );

    label faceI = findMax(magRadSqr);

    if (debug)
    {
        Info<< "findFaceMaxRadius(const pointField&)" << nl
            << "    rotFace  = " << faceI << nl
            << "    point    =  " << faceCentres[faceI] << nl
            << "    distance = " << Foam::sqrt(magRadSqr[faceI])
            << endl;
    }

    return faceI;
}


// * * * * * * * * * * * Protecetd Member Functions  * * * * * * * * * * * * //

void Foam::cyclicAMIPolyPatch::calcTransforms()
{
    if (size())
    {
        // Half0
        const cyclicAMIPolyPatch& half0 = *this;
        vectorField half0Areas(half0.size());
        forAll(half0, facei)
        {
            half0Areas[facei] = half0[facei].normal(half0.points());
        }

        // Half1
        const cyclicAMIPolyPatch& half1 = nbrPatch();
        vectorField half1Areas(half1.size());
        forAll(half1, facei)
        {
            half1Areas[facei] = half1[facei].normal(half1.points());
        }

        calcTransforms
        (
            half0,
            half0.faceCentres(),
            half0Areas,
            half1.faceCentres(),
            half1Areas
        );
    }
}


void Foam::cyclicAMIPolyPatch::calcTransforms
(
    const primitivePatch& half0,
    const pointField& half0Ctrs,
    const vectorField& half0Areas,
    const pointField& half1Ctrs,
    const vectorField& half1Areas
)
{
    if (transform_ != nbrPatch().transform_)
    {
        FatalErrorIn("cyclicAMIPolyPatch::calcTransforms()")
            << "Patch " << name()
            << " has transform type " << transformTypeNames[transform_]
            << ", neighbour patch " << nbrPatchName_
            << " has transform type "
            << nbrPatch().transformTypeNames[nbrPatch().transform_]
            << exit(FatalError);
    }


    // Calculate transformation tensors

    if ((half0Ctrs.size() <= 0) || (half1Ctrs.size() <= 0))
    {
        return;
    }

    switch (transform_)
    {
        case ROTATIONAL:
        {
            label face0 = findFaceMaxRadius(half0Ctrs);
            label face1 = findFaceMaxRadius(half1Ctrs);

            vector n0 = ((half0Ctrs[face0] - rotationCentre_) ^ rotationAxis_);
            vector n1 = ((half1Ctrs[face1] - rotationCentre_) ^ -rotationAxis_);
            n0 /= mag(n0) + VSMALL;
            n1 /= mag(n1) + VSMALL;

            if (debug)
            {
                Pout<< "cyclicAMIPolyPatch::calcTransforms :"
                    << " Specified rotation :"
                    << " n0:" << n0 << " n1:" << n1 << endl;
            }

            // Extended tensor from two local coordinate systems calculated
            // using normal and rotation axis
            const tensor E0
            (
                rotationAxis_,
                (n0 ^ rotationAxis_),
                n0
            );
            const tensor E1
            (
                rotationAxis_,
                (-n1 ^ rotationAxis_),
                -n1
            );
            const tensor revT(E1.T() & E0);
            const_cast<tensorField&>(forwardT()) = tensorField(1, revT.T());
            const_cast<tensorField&>(reverseT()) = tensorField(1, revT);
            const_cast<vectorField&>(separation()).setSize(0);
            const_cast<boolList&>(collocated()) = boolList(1, false);

            break;
        }
        case TRANSLATIONAL:
        {
            // Transform 0 points
            if (debug)
            {
                Pout<< "cyclicAMIPolyPatch::calcTransforms :"
                    << "Specified translation : " << separationVector_
                    << endl;
            }

            const_cast<tensorField&>(forwardT()).clear();
            const_cast<tensorField&>(reverseT()).clear();
            const_cast<vectorField&>(separation()) = vectorField
            (
                1,
                separationVector_
            );
            const_cast<boolList&>(collocated()) = boolList(1, false);

            break;
        }
        default:
        {
            if (debug)
            {
                Pout<< "Assuming cyclic AMI pairs are colocated" << endl;
            }

            const_cast<tensorField&>(forwardT()).clear();
            const_cast<tensorField&>(reverseT()).clear();
            const_cast<vectorField&>(separation()) = vectorField
            (
                1,
                vector::zero
            );
            const_cast<boolList&>(collocated()) = boolList(1, true);

            break;
        }
    }

    if (debug)
    {
        Pout<< "patch: " << name() << nl
            << "    forwardT = " << forwardT() << nl
            << "    reverseT = " << reverseT() << nl
            << "    separation = " << separation() << nl
            << "    collocated = " << collocated() << nl << endl;
    }
}


const Foam::searchableSurface& Foam::cyclicAMIPolyPatch::surf()
{
    word surfType(surfDict_.lookup("type"));
    word surfName(surfDict_.lookupOrDefault("name", name()));

    const polyMesh& mesh = boundaryMesh().mesh();

    surfPtr_ =
        searchableSurface::New
        (
            surfType,
            IOobject
            (
                surfName,
                boundaryMesh().mesh().time().constant(),
                "triSurface",
                boundaryMesh().mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            surfDict_
        );

    return surfPtr_();
}


const Foam::AMIPatchToPatchInterpolation& Foam::cyclicAMIPolyPatch::AMI()
{
    if (!AMIPtr_.valid() && owner())
    {
        const polyPatch& nbr = nbrPatch();
        pointField nbrPoints = nbrPatch().localPoints();

        if (debug)
        {
            OFstream os(name() + "_neighbourPatch-org.obj");
            meshTools::writeOBJ(os, nbrPatch().localFaces(), nbrPoints);
        }

        // transform neighbour patch to local system
        transformPosition(nbrPoints);
        primitivePatch nbrPatch0
        (
            SubList<face>
            (
                nbr.localFaces(),
                nbr.size()
            ),
            nbrPoints
        );

        if (debug)
        {
            OFstream osN(name() + "_neighbourPatch-trans.obj");
            meshTools::writeOBJ(osN, nbrPatch0.localFaces(), nbrPoints);

            OFstream osO(name() + "_ownerPatch.obj");
            meshTools::writeOBJ(osO, this->localFaces(), localPoints());
        }

        bool projectPoints(readBool(surfDict_.lookup("projectPoints")));

        // Construct/apply AMI interpolation to determine addressing and weights
        AMIPtr_.reset
        (
            new AMIPatchToPatchInterpolation
            (
                *this,
                nbrPatch0,
                surf(),
                faceAreaIntersect::tmMesh,
                projectPoints
            )
        );
    }

    return AMIPtr_();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cyclicAMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    polyPatch::initGeometry(pBufs);
}


void Foam::cyclicAMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    calcGeometry
    (
        *this,
        faceCentres(),
        faceAreas(),
        faceCellCentres(),
        nbrPatch().faceCentres(),
        nbrPatch().faceAreas(),
        nbrPatch().faceCellCentres()
    );
}


void Foam::cyclicAMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::initMovePoints(pBufs, p);
}


void Foam::cyclicAMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);

    calcTransforms();

    if (owner())
    {
        AMIPtr_.clear();
        AMI();
    }
}


void Foam::cyclicAMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);
}


void Foam::cyclicAMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    transform_(UNKNOWN),
    rotationAxis_(vector::zero),
    rotationCentre_(point::zero),
    separationVector_(vector::zero),
    AMIPtr_(NULL),
    surfPtr_(NULL),
    surfDict_(dictionary::null)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    nbrPatchName_(dict.lookup("neighbourPatch")),
    nbrPatchID_(-1),
    transform_(UNKNOWN),
    rotationAxis_(vector::zero),
    rotationCentre_(point::zero),
    separationVector_(vector::zero),
    AMIPtr_(NULL),
    surfPtr_(NULL),
    surfDict_(dict.subDict("surface"))
{
    if (nbrPatchName_ == name)
    {
        FatalIOErrorIn
        (
            "cyclicAMIPolyPatch::cyclicAMIPolyPatch"
            "("
                "const word&, "
                "const dictionary&, "
                "const label, "
                "const polyBoundaryMesh&"
            ")",
            dict
        )   << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    if (dict.found("transform"))
    {
        transform_ = transformTypeNames.read(dict.lookup("transform"));
        switch (transform_)
        {
            case ROTATIONAL:
            {
                dict.lookup("rotationAxis") >> rotationAxis_;
                dict.lookup("rotationCentre") >> rotationCentre_;

                scalar magRot = mag(rotationAxis_);
                if (magRot < SMALL)
                {
                    FatalIOErrorIn
                    (
                        "cyclicAMIPolyPatch::cyclicAMIPolyPatch"
                        "("
                            "const word&, "
                            "const dictionary&, "
                            "const label, "
                            "const polyBoundaryMesh&"
                        ")",
                        dict
                    )   << "Illegal rotationAxis " << rotationAxis_ << endl
                        << "Please supply a non-zero vector."
                        << exit(FatalIOError);
                }
                rotationAxis_ /= magRot;

                break;
            }
            case TRANSLATIONAL:
            {
                dict.lookup("separationVector") >> separationVector_;
                break;
            }
            default:
            {
                // no additional info required
            }
        }
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    nbrPatchName_(pp.nbrPatchName_),
    nbrPatchID_(-1),
    transform_(UNKNOWN),
    rotationAxis_(vector::zero),
    rotationCentre_(point::zero),
    separationVector_(vector::zero),
    AMIPtr_(NULL),
    surfPtr_(NULL),
    surfDict_(dictionary::null)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    nbrPatchName_(nbrPatchName),
    nbrPatchID_(-1),
    transform_(UNKNOWN),
    rotationAxis_(vector::zero),
    rotationCentre_(point::zero),
    separationVector_(vector::zero),
    AMIPtr_(NULL),
    surfPtr_(NULL),
    surfDict_(dictionary::null)
{
    if (nbrPatchName_ == name())
    {
        FatalErrorIn
        (
            "const cyclicAMIPolyPatch& "
            "const polyBoundaryMesh&, "
            "const label, "
            "const label, "
            "const label, "
            "const word&"
        )   << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    nbrPatchName_(pp.nbrPatchName_),
    nbrPatchID_(-1),
    transform_(pp.transform_),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    separationVector_(pp.separationVector_),
    AMIPtr_(NULL),
    surfPtr_(NULL),
    surfDict_(pp.surfDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicAMIPolyPatch::~cyclicAMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cyclicAMIPolyPatch::nbrPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        nbrPatchID_ = this->boundaryMesh().findPatchID(nbrPatchName_);

        if (nbrPatchID_ == -1)
        {
            FatalErrorIn("cyclicPolyAMIPatch::nbrPatchID() const")
                << "Illegal neighbourPatch name " << nbrPatchName_
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a cyclic AMI patch
        const cyclicAMIPolyPatch& nbrPatch =
            refCast<const cyclicAMIPolyPatch>
            (
                this->boundaryMesh()[nbrPatchID_]
            );

        if (nbrPatch.nbrPatchName() != name())
        {
            WarningIn("cyclicAMIPolyPatch::nbrPatchID() const")
                << "Patch " << name()
                << " specifies neighbour patch " << nbrPatchName()
                << nl << " but that in return specifies "
                << nbrPatch.nbrPatchName() << endl;
        }
    }

    return nbrPatchID_;
}


bool Foam::cyclicAMIPolyPatch::owner() const
{
    return index() < nbrPatchID();
}


void Foam::cyclicAMIPolyPatch::transformPosition(pointField& l) const
{
    if (!parallel())
    {
        if (transform_ == ROTATIONAL)
        {
            l = Foam::transform(forwardT(), l - rotationCentre_)
              + rotationCentre_;
        }
        else
        {
            l = Foam::transform(forwardT(), l);
        }
    }
    else if (separated())
    {
        // transformPosition gets called on the receiving side,
        // separation gets calculated on the sending side so subtract

        const vectorField& s = separation();
        if (s.size() == 1)
        {
            forAll(l, i)
            {
                l[i] -= s[0];
            }
        }
        else
        {
            l -= s;
        }
    }
}


void Foam::cyclicAMIPolyPatch::transformPosition
(
    point& l,
    const label faceI
) const
{
    if (!parallel())
    {
        const tensor& T =
        (
            forwardT().size() == 1
          ? forwardT()[0]
          : forwardT()[faceI]
        );

        if (transform_ == ROTATIONAL)
        {
            l = Foam::transform(T, l - rotationCentre_) + rotationCentre_;
        }
        else
        {
            l = Foam::transform(T, l);
        }
    }
    else if (separated())
    {
        const vector& s =
        (
            separation().size() == 1
          ? separation()[0]
          : separation()[faceI]
        );

        l -= s;
    }
}


void Foam::cyclicAMIPolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const pointField& thisCtrs,
    const vectorField& thisAreas,
    const pointField& thisCc,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const pointField& nbrCc
)
{
    calcTransforms
    (
        referPatch,
        thisCtrs,
        thisAreas,
        nbrCtrs,
        nbrAreas
    );

    if (owner())
    {
        AMIPtr_.clear();
        AMI();
    }
}


void Foam::cyclicAMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{}


bool Foam::cyclicAMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    // do nothing
    return false;
}


void Foam::cyclicAMIPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);
    os.writeKeyword("neighbourPatch") << nbrPatchName_
        << token::END_STATEMENT << nl;
    switch (transform_)
    {
        case ROTATIONAL:
        {
            os.writeKeyword("transform") << transformTypeNames[transform_]
                << token::END_STATEMENT << nl;
            os.writeKeyword("rotationAxis") << rotationAxis_
                << token::END_STATEMENT << nl;
            os.writeKeyword("rotationCentre") << rotationCentre_
                << token::END_STATEMENT << nl;
            break;
        }
        case TRANSLATIONAL:
        {
            os.writeKeyword("transform") << transformTypeNames[transform_]
                << token::END_STATEMENT << nl;
            os.writeKeyword("separationVector") << separationVector_
                << token::END_STATEMENT << nl;
            break;
        }
        case NOORDERING:
        {
            os.writeKeyword("transform") << transformTypeNames[transform_]
                << token::END_STATEMENT << nl;
            break;
        }
        default:
        {
            // no additional info to write
        }
    }

    os.writeKeyword(surfDict_.dictName());
    os  << surfDict_;
}


// ************************************************************************* //
