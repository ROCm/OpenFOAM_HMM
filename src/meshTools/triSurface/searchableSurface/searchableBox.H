/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::searchableBox

Description
    Searching on bounding box

SourceFiles
    searchableBox.C

\*---------------------------------------------------------------------------*/

#ifndef searchableBox_H
#define searchableBox_H

#include "searchableSurface.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class searchableBox Declaration
\*---------------------------------------------------------------------------*/

class searchableBox
:
    public searchableSurface,
    public treeBoundBox
{
private:

    // Private member data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        searchableBox(const searchableBox&);

        //- Disallow default bitwise assignment
        void operator=(const searchableBox&);


public:

    //- Runtime type information
    TypeName("searchableBox");


    // Constructors

        //- Construct from components
        searchableBox(const word& name, const treeBoundBox& bb);

        searchableBox
        (
            const word& name,
            const objectRegistry& obj,
            const dictionary& dict
        );

    // Destructor

        virtual ~searchableBox();


    // Member Functions

        //- Calculate nearest point on surface. Returns
        //  - bool : any point found nearer than nearestDistSqr
        //  - label: relevant index in surface
        //  - point: actual nearest point found
        virtual pointIndexHit findNearest
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const;

        //- Calculate nearest point on edge. Returns
        //  - bool : any point found nearer than nearestDistSqr
        //  - label: relevant index in surface
        //  - point: actual nearest point found
        virtual pointIndexHit findNearestOnEdge
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const;

        //- Find nearest to line. Returns
        //  - bool : any point found?
        //  - label: relevant index in shapes
        //  - point: actual nearest point found
        //  sets:
        //  - tightest  : bounding box
        //  - linePoint : corresponding nearest point on line
        virtual pointIndexHit findNearest
        (
            const linePointRef& ln,
            treeBoundBox& tightest,
            point& linePoint
        ) const;

        //- Find nearest intersection of line between start and end.
        virtual pointIndexHit findLine
        (
            const point& start,
            const point& end
        ) const;

        //- Find any intersection of line between start and end.
        virtual pointIndexHit findLineAny
        (
            const point& start,
            const point& end
        ) const;

        //- Determine type (inside/outside/mixed) for point. unknown if
        //  cannot be determined (e.g. non-manifold surface)
        virtual volumeType getVolumeType(const point&) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
