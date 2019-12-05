/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::treeDataPoint

Description
    Holds (reference to) pointField. Encapsulation of data needed for
    octree searches.
    Used for searching for nearest point. No bounding boxes around points.
    Only overlaps and calcNearest are implemented, rest makes little sense.

    Optionally works on subset of points.

SourceFiles
    treeDataPoint.C

\*---------------------------------------------------------------------------*/

#ifndef treeDataPoint_H
#define treeDataPoint_H

#include "pointField.H"
#include "treeBoundBox.H"
#include "linePointRef.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
template<class Type> class indexedOctree;

/*---------------------------------------------------------------------------*\
                        Class treeDataPoint Declaration
\*---------------------------------------------------------------------------*/

class treeDataPoint
{
    // Private Data

        const pointField& points_;

        //- Subset of points to work on (or empty)
        const labelList pointLabels_;

        const bool useSubset_;

public:


    class findNearestOp
    {
        const indexedOctree<treeDataPoint>& tree_;

    public:

        findNearestOp(const indexedOctree<treeDataPoint>& tree);

        void operator()
        (
            const labelUList& indices,
            const point& sample,

            scalar& nearestDistSqr,
            label& minIndex,
            point& nearestPoint
        ) const;

        void operator()
        (
            const labelUList& indices,
            const linePointRef& ln,

            treeBoundBox& tightest,
            label& minIndex,
            point& linePoint,
            point& nearestPoint
        ) const;
    };


    class findIntersectOp
    {
    public:

        findIntersectOp(const indexedOctree<treeDataPoint>& tree);

        //- Calculate intersection of triangle with ray.
        //  Sets result accordingly
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    // Declare name of the class and its debug switch
    ClassName("treeDataPoint");


    // Constructors

        //- Construct from pointField
        //  \note Holds reference to the points!
        explicit treeDataPoint(const pointField& points);

        //- Construct from subset of pointField, copies point ids
        //  \note Holds reference to the points!
        treeDataPoint
        (
            const pointField& points,
            const labelUList& pointLabels,
            const bool useSubsetPoints = true
        );

        //- Construct from subset of pointField, moves point ids
        //  \note Holds reference to the points!
        treeDataPoint
        (
            const pointField& points,
            labelList&& pointLabels,
            const bool useSubsetPoints = true
        );


    // Member Functions

        // Access

            //- An empty effective point field?
            inline bool empty() const
            {
                return
                (
                    useSubset_
                  ? pointLabels_.empty()
                  : points_.empty()
                );
            }

            //- The effective point field size
            inline label size() const
            {
                return
                (
                    useSubset_
                  ? pointLabels_.size()
                  : points_.size()
                );
            }

            //- The original point field
            inline const pointField& points() const
            {
                return points_;
            }

            //- The original point ids
            inline const labelList& pointLabels() const
            {
                return pointLabels_;
            }

            //- Use a subset of points
            inline bool useSubset() const
            {
                return useSubset_;
            }

            //- The original (non-subset) point label
            inline label pointLabel(const label index) const
            {
                return
                (
                    useSubset_ && index >= 0
                  ? pointLabels_[index]
                  : index
                );
            }

            //- Point at specified index
            inline const point& shapePoint(const label index) const
            {
                return
                (
                    useSubset_
                  ? points_[pointLabels_[index]]
                  : points_[index]
                );
            }

            //- Representative point cloud for all shapes inside
            //- (one point per shape)
            pointField shapePoints() const;


        // Search

            //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
            //  Only makes sense for closed surfaces.
            volumeType getVolumeType
            (
                const indexedOctree<treeDataPoint>& os,
                const point& sample
            ) const;

            //- Does (bb of) shape at index overlap bb
            bool overlaps
            (
                const label index,
                const treeBoundBox& sampleBb
            ) const;

            //- Does shape at index overlap the sphere
            bool overlaps
            (
                const label index,
                const point& centre,
                const scalar radiusSqr
            ) const;


    // Member Operators

        //- The point at the specified index
        inline const point& operator[](const label index) const
        {
            return shapePoint(index);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
