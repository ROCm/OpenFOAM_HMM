/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "ReferredCell.H"
#include "InteractionLists.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
void Foam::ReferredCell<ParticleType>::setConstructionData
(
    const polyMesh& mesh,
    const label sourceCell
)
{
    // Points

    const labelList& points = mesh.cellPoints()[sourceCell];

    vectorList sourceCellVertices(points.size());

    forAll(sourceCellVertices, sCV)
    {
        sourceCellVertices[sCV] = mesh.points()[points[sCV]];
    }

    vertexPositions_ = referPositions(sourceCellVertices);


    // Edges

    const labelList& edges = mesh.cellEdges()[sourceCell];

    edgeList sourceCellEdges(edges.size());

    forAll(sourceCellEdges, sCE)
    {
        sourceCellEdges[sCE] = mesh.edges()[edges[sCE]];
    }

    locallyMapEdgeList(points, sourceCellEdges);


    // Faces

    labelList faces(mesh.cells()[sourceCell]);

    vectorList sourceCellFaceCentres(faces.size());

    vectorList sourceCellFaceAreas(faces.size());

    labelListList sourceCellFaces(faces.size());

    forAll(faces, f)
    {
        sourceCellFaces[f] = mesh.faces()[faces[f]];

        sourceCellFaceCentres[f] = mesh.faceCentres()[faces[f]];

        sourceCellFaceAreas[f] = mesh.faceAreas()[faces[f]];
    }

    locallyMapFaceList(points, sourceCellFaces);

    faceCentres_ = referPositions(sourceCellFaceCentres);

    faceAreas_ = rotateVectors(sourceCellFaceAreas);
}


template<class ParticleType>
void Foam::ReferredCell<ParticleType>::locallyMapEdgeList
(
    const labelList& points,
    const edgeList& sourceCellEdges
)
{
    edges_.setSize(sourceCellEdges.size());

    forAll(sourceCellEdges, sCE)
    {
        const edge& e(sourceCellEdges[sCE]);

        edges_[sCE].start() = findIndex(points, e.start());

        edges_[sCE].end() = findIndex(points, e.end());

        if
        (
            edges_[sCE].start() == -1
         || edges_[sCE].end() == -1
        )
        {
            FatalErrorIn("Foam::ReferredCell::locallyMapEdgeList")
                << "edgeList and points labelList for "
                << "referred cell do not match: "
                << nl << "points: " << points
                << nl << "egdes: " << sourceCellEdges
                << abort(FatalError);
        }
    }
}


template<class ParticleType>
void Foam::ReferredCell<ParticleType>::locallyMapFaceList
(
    const labelList& points,
    const labelListList& sourceCellFaces
)
{
    faces_.setSize(sourceCellFaces.size());

    forAll(sourceCellFaces, sCF)
    {
        const labelList& sourceCellFace(sourceCellFaces[sCF]);

        labelList& localFace(faces_[sCF]);

        localFace.setSize(sourceCellFace.size());

        forAll(sourceCellFace, p)
        {
            localFace[p] = findIndex(points, sourceCellFace[p]);

            if (localFace[p] == -1)
            {
                FatalErrorIn("Foam::ReferredCell::locallyMapEdgeList")
                    << "edgeList and points labelList for "
                    << "referred cell do not match: "
                    << nl << "points: " << points
                    << nl << "faces: " << sourceCellFaces
                    << abort(FatalError);
            }
        }
    }
}


template<class ParticleType>
Foam::vector Foam::ReferredCell<ParticleType>::referPosition
(
    const vector& positionToRefer
)
{
    return offset_ + (rotation_ & positionToRefer);
}


template<class ParticleType>
Foam::vectorList
Foam::ReferredCell<ParticleType>::referPositions
(
    const vectorList& positionsToRefer
)
{
    return offset_ + (rotation_ & positionsToRefer);
}


template<class ParticleType>
Foam::vector
Foam::ReferredCell<ParticleType>::rotateVector(const vector& vectorToRotate)
{
    return rotation_ & vectorToRotate;
}


template<class ParticleType>
Foam::vectorList
Foam::ReferredCell<ParticleType>::rotateVectors
(
    const vectorList& vectorsToRotate
)
{
    return rotation_ & vectorsToRotate;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::ReferredCell<ParticleType>::ReferredCell()
:
    IDLList<ParticleType>(),
    sourceProc_(-1),
    sourceCell_(-1),
    vertexPositions_(),
    offset_(vector::zero),
    rotation_(I)
{}



template<class ParticleType>
Foam::ReferredCell<ParticleType>::ReferredCell
(
    const polyMesh& mesh,
    const label sourceProc,
    const label sourceCell,
    const vector& offset,
    const tensor& rotation
)
:
    IDLList<ParticleType>(),
    sourceProc_(sourceProc),
    sourceCell_(sourceCell),
    offset_(offset),
    rotation_(rotation)
{
    setConstructionData(mesh, sourceCell);
}


template<class ParticleType>
Foam::ReferredCell<ParticleType>::ReferredCell
(
    const label sourceProc,
    const label sourceCell,
    const vectorList& vertexPositions,
    const edgeList& localEdges,
    const labelListList& localFaces,
    const vectorList& faceCentres,
    const vectorList& faceAreas,
    const vector& offset,
    const tensor& rotation
)
:
    IDLList<ParticleType>(),
    sourceProc_(sourceProc),
    sourceCell_(sourceCell),
    edges_(localEdges),
    faces_(localFaces),
    offset_(offset),
    rotation_(rotation)
{
    // Supplied vertexPositions, faceCentres and faceAreas are of the
    // "original" cell, and need to be transformed to the referred
    // locations on construction

    vertexPositions_ = referPositions(vertexPositions);

    faceCentres_ = referPositions(faceCentres);

    faceAreas_ = rotateVectors(faceAreas);
}


template<class ParticleType>
Foam::ReferredCell<ParticleType>::ReferredCell
(
    const polyMesh& mesh,
    const label sourceProc,
    const label sourceCell,
    const vector& cS,
    const vector& cD,
    const vector& nS,
    const vector& nD
)
:
    IDLList<ParticleType>(),
    sourceProc_(sourceProc),
    sourceCell_(sourceCell)
{
    // It is assumed that the vectors originating from the faces being referred
    // here are correct periodic faces - i.e. they have the same area etc.

    vector nA = -nS/mag(nS);
    vector nB = nD/mag(nD);

    rotation_ = rotationTensor(nA, nB);

    offset_ = cD - (rotation_ & cS);

    // Allow sourceCell = -1 to create a dummy ReferredCell
    // to obtain the transformation

    if(sourceCell >= 0)
    {
        setConstructionData(mesh, sourceCell);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::ReferredCell<ParticleType>::~ReferredCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
Foam::ReferredCell<ParticleType> Foam::ReferredCell<ParticleType>::reRefer
(
    const vector& cS,
    const vector& cD,
    const vector& nS,
    const vector& nD
)
{
    vector nA = -nS/mag(nS);
    vector nB = nD/mag(nD);

    tensor newRotation = rotationTensor(nA, nB);

    vector newOffset = cD - (newRotation & cS);

    tensor reReferredRotation = newRotation & rotation_;

    vector reReferredOffset = newOffset + (newRotation & offset_);

    return ReferredCell
    (
        sourceProc_,
        sourceCell_,
        rotation_.T() & (vertexPositions_ - offset_),
        edges_,
        faces_,
        rotation_.T() & (faceCentres_ - offset_),
        rotation_.T() & (faceAreas_),
        reReferredOffset,
        reReferredRotation
    );
}


template<class ParticleType>
Foam::vector Foam::ReferredCell<ParticleType>::referPosition
(
    const vector& positionToRefer
) const
{
    return offset_ + (rotation_ & positionToRefer);
}


template<class ParticleType>
Foam::vectorList Foam::ReferredCell<ParticleType>::referPosition
(
    const vectorList& positionsToRefer
) const
{
    return offset_ + (rotation_ & positionsToRefer);
}


template<class ParticleType>
Foam::vector Foam::ReferredCell<ParticleType>::rotateVector
(
    const vector& vectorToRotate
) const
{
    return rotation_ & vectorToRotate;
}


template<class ParticleType>
Foam::vectorList Foam::ReferredCell<ParticleType>::rotateVectors
(
    const vectorList& vectorsToRotate
) const
{
    return rotation_ & vectorsToRotate;
}


template<class ParticleType>
void Foam::ReferredCell<ParticleType>::referInParticle
(
    ParticleType* incomingParticlePtr
)
{
    ParticleType& p = *incomingParticlePtr;

    p.position() = referPosition
    (
        p.position()
    );

    p.transformProperties(rotation_);

    this->append(incomingParticlePtr);
}


template<class ParticleType>
bool Foam::ReferredCell<ParticleType>::duplicate
(
    const ReferredCell<ParticleType>& refCellDupl
) const
{
    return
    (
        sourceProc_ == refCellDupl.sourceProc()
     && sourceCell_ == refCellDupl.sourceCell()
     && mag(offset_ - refCellDupl.offset())
      < InteractionLists<ParticleType>::transTol
    );
}


template<class ParticleType>
bool Foam::ReferredCell<ParticleType>::duplicate
(
    const label procNo,
    const label nCells
) const
{
    return
    (
        sourceProc_ == procNo
     && sourceCell_ < nCells
     && mag(offset_)
      < InteractionLists<ParticleType>::transTol
    );
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class ParticleType>
bool Foam::operator==
(
    const ReferredCell<ParticleType>& a,
    const ReferredCell<ParticleType>& b
)
{
    return const_cast<ReferredCell<ParticleType>&>(a).duplicate
    (
        const_cast<const ReferredCell<ParticleType>&>(b)
    );
}


template<class ParticleType>
bool Foam::operator!=
(
    const ReferredCell<ParticleType>& a,
    const ReferredCell<ParticleType>& b
)
{
    return !(a == b);
}


template<class ParticleType>
Foam::Istream& Foam::operator>>(Istream& is, ReferredCell<ParticleType>& rC)
{

    is  >> rC.sourceProc_
        >> rC.sourceCell_
        >> rC.vertexPositions_
        >> rC.edges_
        >> rC.faces_
        >> rC.faceCentres_
        >> rC.faceAreas_
        >> rC.offset_
        >> rC.rotation_;

    is.check
    (
        "Istream& operator<<(Istream& f, const ReferredCell<ParticleType>& rC"
    );

    return is;
}


template<class ParticleType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReferredCell<ParticleType>& rC
)
{

    os  << rC.sourceProc()
        << token::SPACE << rC.sourceCell()
        << token::SPACE << rC.vertexPositions()
        << token::SPACE << rC.edges()
        << token::SPACE << rC.faces()
        << token::SPACE << rC.faceCentres()
        << token::SPACE << rC.faceAreas()
        << token::SPACE << rC.offset()
        << token::SPACE << rC.rotation();

    os.check
    (
        "Ostream& operator<<(Ostream& f, const ReferredCell<ParticleType>& rC"
    );

    return os;
}


// ************************************************************************* //
