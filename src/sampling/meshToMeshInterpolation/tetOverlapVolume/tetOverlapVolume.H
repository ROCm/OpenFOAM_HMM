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


Class
    Foam::tetOverlapVolume

Description
    Calculates overlap volume of two tets.

SourceFiles
    tetOverlapVolume.C

\*---------------------------------------------------------------------------*/

#ifndef tetOverlapVolume_H
#define tetOverlapVolume_H

#include "tetrahedron.H"
#include "fvMesh.H"
#include "plane.H"
#include "tetPointRef.H"
#include "OFstream.H"
#include "meshTools.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "tetPoints.H"
#include "tetCell.H"
#include "EdgeMap.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class tetOverlapVolume Declaration
\*---------------------------------------------------------------------------*/

class tetOverlapVolume
{
    // Private member functions

        //- Plane intersection
        inline point planeIntersection
        (
            const FixedList<scalar, 4>& d,
            const tetPoints& t,
            const label negI,
            const label posI
        );


        //- Decompose prism
        template <class TetOp> inline void decomposePrism
        (
            const FixedList<point, 6>& points,
            TetOp& op
        );


        //- Helping cľasses
        class dummyTetOp
        {
         public:

            inline void operator()(const tetPoints&){}
        };


        class sumTetVolOp
        {
         public:
            scalar vol_;

            inline sumTetVolOp()
            :
                vol_(0.0)
            {}

            inline void operator()(const tetPoints& tet)
            {
                vol_ += tet.tet().mag();
            }
        };


        class storeTetOp
        {
            FixedList<tetPoints, 200>& tets_;
            label& nTets_;

         public:

            inline storeTetOp(FixedList<tetPoints, 200>& tets, label& nTets)
            :
                tets_(tets),
                nTets_(nTets)
            {}

            inline void operator()(const tetPoints& tet)
            {
                tets_[nTets_++] = tet;
            }
        };


        //- Slice. Split tet into subtets above and below plane
        template <class AboveTetOp, class BelowTetOp>
        inline void tetSliceWithPlane
        (
            const tetPoints& tet,
            const plane& pl,

            AboveTetOp& aboveOp,
            BelowTetOp& belowOp
        );


        //- Tet overlap
        void tetTetOverlap
        (
            const tetPoints& tetA,
            const tetPoints& tetB,
            FixedList<tetPoints, 200>& insideTets,
            label& nInside,
            FixedList<tetPoints, 200>& outsideTets,
            label& nOutside
        );


        //- Tet Overlap Vol
        inline scalar tetTetOverlapVol
        (
            const tetPoints& tetA,
            const tetPoints& tetB
        );


        //- Return a const treeBoundBox
        inline const treeBoundBox pyrBb
        (
            const pointField& points,
            const face& f,
            const point& fc
        );


public:


    //- Runtime type information
    TypeName("tetOverlapVolume");


    // Constructors

        //- Null constructor
        tetOverlapVolume();


    //- Destructor
     virtual ~tetOverlapVolume();


    // Public members

        //- Return a list of cells in meshA which overlaps with cellBI in
        // meshB
        labelList overlappingCells
        (
            const fvMesh& meshA,
            const fvMesh& meshB,
            const label cellBI
        ) const;


        //- Calculates the overlap volume
        scalar cellCellOverlapVolumeMinDecomp
        (
            const primitiveMesh& meshA,
            const label cellAI,

            const primitiveMesh& meshB,
            const label cellBI,
            const treeBoundBox& cellBbB
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
