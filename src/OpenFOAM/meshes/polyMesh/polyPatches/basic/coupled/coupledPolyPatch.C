/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "coupledPolyPatch.H"
#include "SortableList.H"
#include "ListOps.H"
#include "transform.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledPolyPatch, 0);
}

Foam::scalar Foam::coupledPolyPatch::matchTol_ = 1E-3;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::coupledPolyPatch::writeOBJ(Ostream& os, const point& pt)
{
    os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
}


void Foam::coupledPolyPatch::writeOBJ
(
    Ostream& os,
    const pointField& points,
    const labelList& pointLabels
)
{
    forAll(pointLabels, i)
    {
        writeOBJ(os, points[pointLabels[i]]);
    }
}


void Foam::coupledPolyPatch::writeOBJ
(
    Ostream& os,
    const point& p0,
    const point& p1,
    label& vertI
)
{
    writeOBJ(os, p0);
    vertI++;

    writeOBJ(os, p1);
    vertI++;

    os<< "l " << vertI-1 << ' ' << vertI << nl;
}


void Foam::coupledPolyPatch::writeOBJ
(
    const fileName& fName,
    const faceList& faces,
    const pointField& points
)
{
    OFstream os(fName);

    Map<label> foamToObj(4*faces.size());

    label vertI = 0;

    forAll(faces, i)
    {
        const face& f = faces[i];

        forAll(f, fp)
        {
            if (foamToObj.insert(f[fp], vertI))
            {
                writeOBJ(os, points[f[fp]]);
                vertI++;
            }
        }

        os << 'l';
        forAll(f, fp)
        {
            os << ' ' << foamToObj[f[fp]]+1;
        }
        os << ' ' << foamToObj[f[0]]+1 << nl;
    }
}


Foam::pointField Foam::coupledPolyPatch::calcFaceCentres
(
    const faceList& faces,
    const pointField& points
)
{
    pointField ctrs(faces.size());

    forAll(faces, faceI)
    {
        ctrs[faceI] = faces[faceI].centre(points);
    }

    return ctrs;
}


Foam::pointField Foam::coupledPolyPatch::getAnchorPoints
(
    const faceList& faces,
    const pointField& points
)
{
    pointField anchors(faces.size());

    forAll(faces, faceI)
    {
        anchors[faceI] = points[faces[faceI][0]];
    }

    return anchors;
}


bool Foam::coupledPolyPatch::inPatch
(
    const labelList& oldToNew,
    const label oldFaceI
) const
{
    label faceI = oldToNew[oldFaceI];

    return faceI >= start() && faceI < start()+size();
}


Foam::label Foam::coupledPolyPatch::whichPatch
(
    const labelList& patchStarts,
    const label faceI
)
{
    forAll(patchStarts, patchI)
    {
        if (patchStarts[patchI] <= faceI)
        {
            if (patchI == patchStarts.size()-1)
            {
                return patchI;
            }
            else if (patchStarts[patchI+1] > faceI)
            {
                return patchI;
            }
        }
    }

    return -1;
}


Foam::scalarField Foam::coupledPolyPatch::calcFaceTol
(
    const faceList& faces,
    const pointField& points,
    const pointField& faceCentres
)
{
    // Calculate typical distance per face
    scalarField tols(faces.size());

    forAll(faces, faceI)
    {
        const point& cc = faceCentres[faceI];

        const face& f = faces[faceI];

        scalar maxLen = -GREAT;

        forAll(f, fp)
        {
            maxLen = max(maxLen, mag(points[f[fp]] - cc));
        }
        tols[faceI] = matchTol_ * maxLen;
    }
    return tols;
}


Foam::label Foam::coupledPolyPatch::getRotation
(
    const pointField& points,
    const face& f,
    const point& anchor,
    const scalar tol
)
{
    label anchorFp = -1;
    scalar minDistSqr = GREAT;

    forAll(f, fp)
    {
        scalar distSqr = magSqr(anchor - points[f[fp]]);

        if (distSqr < minDistSqr)
        {
            minDistSqr = distSqr;
            anchorFp = fp;
        }
    }

    if (anchorFp == -1 || mag(minDistSqr) > tol)
    {
        return -1;
    }
    else
    {
        // Positive rotation
        return (f.size() - anchorFp) % f.size();
    }
}


void Foam::coupledPolyPatch::calcTransformTensors
(
    const vector& Cf,
    const vector& Cr,
    const vector& nf,
    const vector& nr
) const
{
    if (debug)
    {
        Pout<< "coupledPolyPatch::calcTransformTensors : " << name() << endl
            << "    nf:" << nf << nl
            << "    nr:" << nr << nl
            << "    mag(nf&nr):" << mag(nf & nr) << nl
            << "    Cf:" << Cf << nl
            << "    Cr:" << Cr << endl;
    }

    if (mag(nf & nr) < 1 - SMALL)
    {
        separation_.setSize(0);

        forwardT_ = tensorField(1, rotationTensor(-nr, nf));
        reverseT_ = tensorField(1, rotationTensor(nf, -nr));
    }
    else
    {
        forwardT_.setSize(0);
        reverseT_.setSize(0);

        vector separation = (nf & (Cr - Cf))*nf;

        if (mag(separation) > SMALL)
        {
            separation_ = vectorField(1, separation);
        }
        else
        {
            separation_.setSize(0);
        }
    }

    if (debug)
    {
        Pout<< "    separation_:" << separation_ << nl
            << "    forwardT size:" << forwardT_.size() << endl;
    }
}


void Foam::coupledPolyPatch::calcTransformTensors
(
    const vectorField& Cf,
    const vectorField& Cr,
    const vectorField& nf,
    const vectorField& nr
) const
{
    if (debug)
    {
        Pout<< "coupledPolyPatch::calcTransformTensors : " << name() << endl
            << "    size:" << size() << nl
            << "    sum(mag(nf & nr)):" << sum(mag(nf & nr)) << endl;
    }

    if (size() == 0)
    {
        // Dummy geometry.
        separation_.setSize(0);
        forwardT_ = I;
        reverseT_ = I;
    }
    else if (sum(mag(nf & nr)) < size() - 1E-12)
    {
        separation_.setSize(0);

        forwardT_.setSize(size());
        reverseT_.setSize(size());

        forAll (forwardT_, facei)
        {
            forwardT_[facei] = rotationTensor(-nr[facei], nf[facei]);
            reverseT_[facei] = rotationTensor(nf[facei], -nr[facei]);
        }

        if (sum(mag(forwardT_ - forwardT_[0])) < 1E-12)
        {
            forwardT_.setSize(1);
            reverseT_.setSize(1);
        }
    }
    else
    {
        forwardT_.setSize(0);
        reverseT_.setSize(0);

        separation_ = (nf&(Cr - Cf))*nf;

        if (sum(mag(separation_))/size() < 1E-12)
        {
            separation_.setSize(0);
        }
        else if (sum(mag(separation_ - separation_[0]))/size() < 1E-12)
        {
            separation_.setSize(1);
        }
    }

    if (debug)
    {
        Pout<< "    separation_:" << separation_ << nl
            << "    forwardT size:" << forwardT_.size() << endl;
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::coupledPolyPatch::coupledPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm)
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm)
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const coupledPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm)
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const coupledPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledPolyPatch::~coupledPolyPatch()
{}


// ************************************************************************* //
