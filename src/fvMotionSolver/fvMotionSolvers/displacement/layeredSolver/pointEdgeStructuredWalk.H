/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::pointEdgeStructuredWalk

Description
    Determines length of string of edges walked to point.

SourceFiles
    pointEdgeStructuredWalkI.H
    pointEdgeStructuredWalk.C

\*---------------------------------------------------------------------------*/

#ifndef pointEdgeStructuredWalk_H
#define pointEdgeStructuredWalk_H

#include "point.H"
#include "label.H"
#include "scalar.H"
#include "tensor.H"
#include "pTraits.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class pointEdgeStructuredWalk Declaration
\*---------------------------------------------------------------------------*/

class pointEdgeStructuredWalk
{
    // Private data

        //- Starting location
        point point0_;

        //- Previous point
        point previousPoint_;

        //- Sum of distance
        scalar dist_;

        //- Passive data
        vector data_;

    // Private Member Functions

        //- Evaluate distance to point.
        inline bool update
        (
            const pointEdgeStructuredWalk& w2,
            const scalar tol
        );

public:

    // Constructors

        //- Construct null
        inline pointEdgeStructuredWalk();

        //- Construct from components
        inline pointEdgeStructuredWalk
        (
            const point&,
            const point&,
            const scalar,
            const vector&
        );


    // Member Functions

        // Access

            inline bool inZone() const;

            inline scalar dist() const;

            inline const vector& data() const;


        // Needed by meshWave

            //- Check whether still contains original (invalid) value.
            inline bool valid() const;

            //- Check for identical geometrical data. Used for cyclics checking.
            inline bool sameGeometry
            (
                const pointEdgeStructuredWalk&,
                const scalar tol
            ) const;

            //- Convert origin to relative vector to leaving point
            //  (= point coordinate)
            inline void leaveDomain
            (
                const polyPatch& patch,
                const label patchPointI,
                const point& pos
            );

            //- Convert relative origin to absolute by adding entering point
            inline void enterDomain
            (
                const polyPatch& patch,
                const label patchPointI,
                const point& pos
            );

            //- Apply rotation matrix to origin
            inline void transform(const tensor& rotTensor);

            //- Influence of edge on point
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointI,
                const label edgeI,
                const pointEdgeStructuredWalk& edgeInfo,
                const scalar tol
            );

            //- Influence of different value on same point.
            //  Merge new and old info.
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointI,
                const pointEdgeStructuredWalk& newPointInfo,
                const scalar tol
            );

            //- Influence of different value on same point.
            //  No information about current position whatsoever.
            inline bool updatePoint
            (
                const pointEdgeStructuredWalk& newPointInfo,
                const scalar tol
            );

            //- Influence of point on edge.
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgeI,
                const label pointI,
                const pointEdgeStructuredWalk& pointInfo,
                const scalar tol
            );


    // Member Operators

        //Note: Used to determine whether to call update.
        inline bool operator==(const pointEdgeStructuredWalk&) const;

        inline bool operator!=(const pointEdgeStructuredWalk&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const pointEdgeStructuredWalk&);
        friend Istream& operator>>(Istream&, pointEdgeStructuredWalk&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointEdgeStructuredWalkI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
