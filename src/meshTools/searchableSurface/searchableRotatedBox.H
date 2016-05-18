/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::searchableRotatedBox

Description
    Searching on a rotated box

    Box defined as min and max coordinate. Rotation by coordinate system
    at box middle.

    E.g. box with sides 1 1 1 rotated 45 degrees around z-axis at
         origin (0.5 0.5 0.5)
    \verbatim
        span     (1 1 1);
        origin  (0.5 0.5 0.5);
        e1      (1 1 0);
        e3      (0 0 1);
    \endverbatim

SourceFiles
    searchableRotatedBox.C

\*---------------------------------------------------------------------------*/

#ifndef searchableRotatedBox_H
#define searchableRotatedBox_H

#include "searchableSurface.H"
#include "coordinateSystem.H"
#include "searchableBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class searchableRotatedBox Declaration
\*---------------------------------------------------------------------------*/

class searchableRotatedBox
:
    public searchableSurface
{
private:

    // Private Member Data

        //- box in local coordinate system
        searchableBox box_;

        //- transformation from local to global coordinates
        coordinateSystem transform_;

        //- (global) corner points (in treeBoundBox order)
        pointField points_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        searchableRotatedBox(const searchableRotatedBox&);

        //- Disallow default bitwise assignment
        void operator=(const searchableRotatedBox&);


public:

    //- Runtime type information
    TypeName("searchableRotatedBox");


    // Constructors

        //- Construct from dictionary (used by searchableSurface)
        searchableRotatedBox
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableRotatedBox();


    // Member Functions

        virtual const wordList& regions() const;

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return true;
        }

        //- Range of local indices that can be returned.
        virtual label size() const
        {
            return 6;
        }

        //- Get representative set of element coordinates
        //  Usually the element centres (should be of length size()).
        virtual tmp<pointField> coordinates() const;

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const;

        // Does any part of the surface overlap the supplied bound box?
        virtual bool overlaps(const boundBox& bb) const;

        // Single point queries.

            //- Inherit findNearest from searchableSurface
            using searchableSurface::findNearest;

            //- Calculate nearest point on surface.
            //  Returns
            //  - bool : any point found nearer than nearestDistSqr
            //  - label: relevant index in surface (=face 0..5)
            //  - point: actual nearest point found
            pointIndexHit findNearest
            (
                const point& sample,
                const scalar nearestDistSqr
            ) const;

            //- Find nearest to segment.
            //  Returns
            //  - bool : any point found?
            //  - label: relevant index in shapes (=face 0..5)
            //  - point: actual nearest point found
            //  sets:
            //  - tightest  : bounding box
            //  - linePoint : corresponding nearest point on line
            pointIndexHit findNearest
            (
                const linePointRef& ln,
                treeBoundBox& tightest,
                point& linePoint
            ) const;

            //- Find nearest intersection of line between start and end.
            pointIndexHit findLine
            (
                const point& start,
                const point& end
            ) const;

            //- Find any intersection of line between start and end.
            pointIndexHit findLineAny
            (
                const point& start,
                const point& end
            ) const;


        // Multiple point queries.

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField&,
                List<volumeType>&
            ) const;


        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                NotImplemented;
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
