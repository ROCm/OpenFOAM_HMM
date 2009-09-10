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

Class
    Foam::ReferredCell

Description

SourceFiles
    ReferredCellI.H
    ReferredCell.C

\*---------------------------------------------------------------------------*/

#ifndef ReferredCell_H
#define ReferredCell_H

#include "vector.H"
#include "vectorList.H"
#include "tensor.H"
#include "transform.H"
#include "IDLList.H"
#include "labelList.H"
#include "edgeList.H"
#include "polyMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class ParticleType>
class ReferredCell;

template<class ParticleType>
bool operator==
(
    const ReferredCell<ParticleType>& a,
    const ReferredCell<ParticleType>& b
);

template<class ParticleType>
bool operator!=
(
    const ReferredCell<ParticleType>& a,
    const ReferredCell<ParticleType>& b
);

template<class ParticleType>
Istream& operator>>
(
    Istream&,
    ReferredCell<ParticleType>&
);

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const ReferredCell<ParticleType>&
);


/*---------------------------------------------------------------------------*\
                          Class ReferredCell Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class ReferredCell
:
    public IDLList<ParticleType>
{
    // Private data

        label sourceProc_;

        label sourceCell_;

        //- Referred vertex positions
        vectorList vertexPositions_;

        edgeList edges_;

        labelListList faces_;

        vectorList faceCentres_;

        vectorList faceAreas_;

        labelList realCellsForInteraction_;

        vector offset_;

        tensor rotation_;


    // Private Member Functions

        void setConstructionData
        (
            const polyMesh& mesh,
            const label sourceCell
        );

        void locallyMapEdgeList
        (
            const labelList& points,
            const edgeList& sourceCellEdges
        );

        void locallyMapFaceList
        (
            const labelList& points,
            const labelListList& sourceCellFaces
        );

        vector referPosition(const vector& positionToRefer);

        vectorList referPositions(const vectorList& positionsToRefer);

        vector rotateVector(const vector& vectorToRotate);

        vectorList rotateVectors(const vectorList& vectorsToRotate);


public:

    // Constructors

        //- Construct null

        ReferredCell();

        //- Construct from components with external edge information
        ReferredCell
        (
            const polyMesh& mesh,
            const label sourceProc,
            const label sourceCell,
            const vector& offset,
            const tensor& rotation
        );

        //- Construct from components with existing local edge information
        ReferredCell
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
        );

        //- Construct from pair of face centers (c) and plain
        //  face normals (n) (no need to make unit vectors or
        //  reverse one direction)
        //  Order of vectors important (S = source, D = Destination).
        //  External edge information.

        ReferredCell
        (
            const polyMesh& mesh,
            const label sourceProc,
            const label sourceCell,
            const vector& cS,
            const vector& cD,
            const vector& nS,
            const vector& nD
        );


    // Destructor

        virtual ~ReferredCell();


    // Member Functions

        //- Take this ReferredCell object that has already had it's transform
        //  calculated and refer it on again, retaining same source info.
        ReferredCell reRefer
        (
            const vector& cS,
            const vector& cD,
            const vector& nS,
            const vector& nD
        );

        //- Use internal transformatation values to transform the given
        //  postion to its new location.
        vector referPosition(const vector& positionToRefer) const;

        //- Use internal transformatation values to transform the given
        //  list of postions to their new locations.
        vectorList referPosition(const vectorList& positionsToRefer) const;

        //- Use internal transformatation values to rotate the given vector
        vector rotateVector(const vector& vectorToRotate) const;

        //- Use internal transformatation values to rotate the given
        // list of vectors
        vectorList rotateVectors(const vectorList& vectorsToRotate) const;

        //- referInParticle takes a referred particle from a source
        //  processor, transforming its position and properties
        void referInParticle(ParticleType* incomingParticlePtr);

        //- duplicate() function to test whether a referred or real cell
        //  supplied by arguement is a duplicate of this ReferredCell.
        //  Can be used bi-directionally - i.e. can be called on an existing
        //  referred cell with a proposed ReferredCell as argument,
        //  or vice versa.  Can only be called by a proposed ReferredCell with
        //  a real cell index as arguement to test to see if the proposed
        //  ReferredCell is a duplicate.
        //  A duplicate cell is defined as one which has the same source
        //  processor, source cell, and an equal offset. Real cells have zero
        //  offset by definition.
        bool duplicate(const ReferredCell& refCellDupl) const;

        bool duplicate(const label procNo, const label nCells) const;


        // Access

            inline label sourceProc() const;

            inline label sourceCell() const;

            inline const vector& offset() const;

            inline const tensor& rotation() const;

            inline const vectorList& vertexPositions() const;

            inline const edgeList& edges() const;

            inline const labelListList& faces() const;

            inline const vectorList& faceCentres() const;

            inline const vectorList& faceAreas() const;

            inline labelList& realCells();

            inline const labelList& realCellsForInteraction() const;


    // Friend Operators

        friend bool operator== <ParticleType>
        (
            const ReferredCell<ParticleType>& a,
            const ReferredCell<ParticleType>& b
        );

        friend bool operator!= <ParticleType>
        (
            const ReferredCell<ParticleType>& a,
            const ReferredCell<ParticleType>& b
        );

        friend Istream& operator>> <ParticleType>
        (
            Istream&,
            ReferredCell<ParticleType>&
        );

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const ReferredCell<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReferredCellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReferredCell.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
