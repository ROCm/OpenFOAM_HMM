/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchPointEdgeCirculator::setEnd()
{
    edgeID_ = -1;
    pointID_ = -1;
}


// Cross face to other edge on point
void Foam::patchPointEdgeCirculator::otherEdge(const label faceI)
{
    const labelList& fEdges = patch_.faceEdges()[faceI];
    const face& f = patch_.localFaces()[faceI];
    label fp = findIndex(f, pointID_);

    if (fEdges[fp] == edgeID_)
    {
        edgeID_ = fEdges[f.rcIndex(fp)];
    }
    else
    {
        // Check for now
        if (fEdges[f.rcIndex(fp)] != edgeID_)
        {
            FatalErrorIn("patchPointEdgeCirculator::otherEdge(const label)")
                << "face:" << faceI
                << " verts:" << f
                << " edges:" << fEdges
                << " looking for other edge than " << edgeID_
                << abort(FatalError);
        }
        edgeID_ = fEdges[fp];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
Foam::patchPointEdgeCirculator::patchPointEdgeCirculator
(
    const primitiveFacePatch& patch,
    const PackedBoolList& nonManifoldEdge,
    const label edgeID,
    const label index,

    const label pointID
)
:
    patch_(patch),
    nonManifoldEdge_(nonManifoldEdge),
    edgeID_(edgeID),
    index_(index),
    pointID_(pointID),
    startEdgeID_(edgeID_)
{
    if (edgeID_ != -1)
    {
        const edge& e = patch_.edges()[edgeID_];

        if (pointID_ != e[0] && pointID_ != e[1])
        {
            FatalErrorIn
            (
                "patchPointEdgeCirculator::patchPointEdgeCirculator(..)"
            )   << "edge " << edgeID_ << " verts " << e
                << " does not contain point " << pointID_ << abort(FatalError);
        }
    }
}


//- Construct copy
Foam::patchPointEdgeCirculator::patchPointEdgeCirculator
(
    const patchPointEdgeCirculator& circ
)
:
    patch_(circ.patch_),
    nonManifoldEdge_(circ.nonManifoldEdge_),
    edgeID_(circ.edgeID_),
    index_(circ.index_),
    pointID_(circ.pointID_),
    startEdgeID_(circ.startEdgeID_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::patchPointEdgeCirculator::edgeID() const
{
    return edgeID_;
}


Foam::label Foam::patchPointEdgeCirculator::index() const
{
    return index_;
}


Foam::label Foam::patchPointEdgeCirculator::pointID() const
{
    return pointID_;
}


Foam::label Foam::patchPointEdgeCirculator::faceID() const
{
    if (edgeID_ != -1 && index_ != -1)
    {
        return patch_.edgeFaces()[edgeID_][index_];
    }
    else
    {
        return -1;
    }
}


void Foam::patchPointEdgeCirculator::operator=
(
    const patchPointEdgeCirculator& circ
)
{
    edgeID_ = circ.edgeID_;
    index_ = circ.index_;
    pointID_ = circ.pointID_;
    startEdgeID_ = circ.startEdgeID_;
}


bool Foam::patchPointEdgeCirculator::operator==
(
    const patchPointEdgeCirculator& circ
) const
{
    // Do just enough to have setEnd() produce an iterator equal to end().
    // Could include the direction(i.e. index_) to make sure that two
    // circulators around same point but in different direction are not equal.
    return edgeID_ == circ.edgeID_ && pointID_ == circ.pointID_;
}


bool Foam::patchPointEdgeCirculator::operator!=
(
    const patchPointEdgeCirculator& circ
) const
{
    return !(*this == circ);
}


void Foam::patchPointEdgeCirculator::setCanonical()
{
    if (edgeID_ == -1)
    {
        FatalErrorIn("patchPointEdgeCirculator::setCanonical()")
            << "Already reached end(). Cannot walk any further."
            << abort(FatalError);
    }

    label minEdgeID = edgeID_;
    label minIndex = index_;

    while (true)
    {
        if (nonManifoldEdge_[edgeID_])
        {
            break;
        }

        // Step back
        const labelList& eFaces = patch_.edgeFaces()[edgeID_];

        if (eFaces.size() != 2)
        {
            FatalErrorIn("patchPointEdgeCirculator::setCanonical()")
                << "problem eFaces:" << eFaces << abort(FatalError);
        }

        label faceI = (index_ == 0 ? eFaces[1] : eFaces[0]);

        // Step to other edge on face
        otherEdge(faceI);

        // Update index
        index_ = findIndex(patch_.edgeFaces()[edgeID_], faceI);

        if (edgeID_ < minEdgeID)
        {
            minEdgeID = edgeID_;
            minIndex = index_;
        }

        if (edgeID_ == startEdgeID_)
        {
            edgeID_ = minEdgeID;
            index_ = minIndex;
            break;
        }
    }

    startEdgeID_ = edgeID_;
}


//- Step to next edge.
Foam::patchPointEdgeCirculator&
Foam::patchPointEdgeCirculator::operator++()
{
    if (index_ == -1)
    {
        setEnd();
    }
    else
    {
        // Step to other edge on face
        label faceI = patch_.edgeFaces()[edgeID_][index_];
        otherEdge(faceI);

        if (edgeID_ == startEdgeID_)
        {
            setEnd();
        }
        else if (nonManifoldEdge_[edgeID_])
        {
            // Reached non-manifold edge. Cannot walk further.
            // Mark so it gets set to end next time.
            index_ = -1;
        }
        else
        {
            const labelList& eFaces = patch_.edgeFaces()[edgeID_];

            if (eFaces.size() != 2)
            {
                FatalErrorIn("patchPointEdgeCirculator:::operator++()")
                    << "problem eFaces:" << eFaces << abort(FatalError);
            }
            // Point to the face that is not faceI
            index_ = (eFaces[0] != faceI ? 0 : 1);
        }
    }

    return *this;
}


Foam::patchPointEdgeCirculator Foam::patchPointEdgeCirculator::begin() const
{
    patchPointEdgeCirculator iter(*this);
    iter.setCanonical();

    return iter;
}


Foam::patchPointEdgeCirculator Foam::patchPointEdgeCirculator::cbegin() const
{
    patchPointEdgeCirculator iter(*this);
    iter.setCanonical();

    return iter;
}


const Foam::patchPointEdgeCirculator& Foam::patchPointEdgeCirculator::end()
 const
{
    return endConstIter;
}

const Foam::patchPointEdgeCirculator& Foam::patchPointEdgeCirculator::cend()
 const
{
    return endConstIter;
}


// ************************************************************************* //
