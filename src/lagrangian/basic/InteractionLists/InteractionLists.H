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
    Foam::InteractionLists

Description

SourceFiles
    InteractionListsI.H
    InteractionLists.C
    InteractionListsIO.C

\*---------------------------------------------------------------------------*/

#ifndef InteractionLists_H
#define InteractionLists_H

#include "polyMesh.H"
#include "vector2D.H"
#include "DirectInteractionList.H"
#include "ReferredCell.H"
#include "ReferredCellList.H"
#include "sendingReferralList.H"
#include "receivingReferralList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class InteractionLists Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class InteractionLists
{
    // Private data

        const polyMesh& mesh_;

        scalar maxDistanceSqr_;

        DirectInteractionList<ParticleType> dil_;

        ReferredCellList<ParticleType> ril_;

        List<sendingReferralList> cellSendingReferralLists_;

        List<receivingReferralList> cellReceivingReferralLists_;


    // Private Member Functions

        //- Build referralLists which define how to send information
        // to ReferredCells to source cells
        void buildCellReferralLists();

        //- Disallow default bitwise copy construct
        InteractionLists(const InteractionLists&);

        //- Disallow default bitwise assignment
        void operator=(const InteractionLists&);


public:

    // Static data members

        //- Tolerance for checking that faces on a patch segment
        static scalar transTol;


    // Constructors

        //- Construct and create all information from the mesh
        InteractionLists
        (
            const polyMesh& mesh,
            scalar maxDistanceSqr,
            bool pointPointListBuild = false
        );

        //- Construct from file
        InteractionLists(const polyMesh& mesh);


    // Destructor

        ~InteractionLists();


    // Member Functions

            bool testPointPointDistance
            (
                const label ptI,
                const label ptJ
            ) const;

            bool testPointFaceDistance
            (
                const label p,
                const label faceNo
            ) const;

            bool testPointFaceDistance
            (
                const label p,
                const ReferredCell<ParticleType>& refCell
            ) const;

            bool testPointFaceDistance
            (
                const vectorList& pointsToTest,
                const label faceNo
            ) const;

            bool testPointFaceDistance
            (
                const vector& p,
                const label faceNo
            ) const;

            bool testPointFaceDistance
            (
                const vector& p,
                const labelList& faceToTest,
                const vectorList& points,
                const vector& faceC,
                const vector& faceA
            ) const;

            bool testEdgeEdgeDistance
            (
                const edge& eI,
                const edge& eJ
            ) const;

            bool testEdgeEdgeDistance
            (
                const edge& eI,
                const vector& eJs,
                const vector& eJe
            ) const;

            const labelList realCellsInRangeOfSegment
            (
                const labelList& segmentFaces,
                const labelList& segmentEdges,
                const labelList& segmentPoints
            ) const;

            const labelList ReferredCellsInRangeOfSegment
            (
                const List<ReferredCell<ParticleType> >&
                    referredInteractionList,
                const labelList& segmentFaces,
                const labelList& segmentEdges,
                const labelList& segmentPoints
            ) const;


        // Access

            inline const polyMesh& mesh() const;

            inline const DirectInteractionList<ParticleType>& dil() const;

            inline const ReferredCellList<ParticleType>& ril() const;
            inline ReferredCellList<ParticleType>& ril();

            inline const List<sendingReferralList>&
                cellSendingReferralLists() const;

            inline const List<receivingReferralList>&
                cellReceivingReferralLists() const;

            inline label nInteractingProcs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InteractionListsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "InteractionLists.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
