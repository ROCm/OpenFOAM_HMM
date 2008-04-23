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

Class
    Foam::searchableSurface

Description
    Base class of (analytical or triangulated) surface.
    Encapsulates all the search routines.

SourceFiles
    searchableSurface.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurface_H
#define searchableSurface_H

#include "pointField.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "pointIndexHit.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class treeBoundBox;

/*---------------------------------------------------------------------------*\
                           Class searchableSurface Declaration
\*---------------------------------------------------------------------------*/

class searchableSurface
{
public:

    // Data types

        //- volume types
        enum volumeType
        {
            UNKNOWN = 0,
            MIXED = 1,
            INSIDE = 2,
            OUTSIDE = 3
        };

private:

    // Private data

        const word name_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        searchableSurface(const searchableSurface&);

        //- Disallow default bitwise assignment
        void operator=(const searchableSurface&);


public:

    //- Runtime type information
    TypeName("searchableSurface");

    // Declare run-time constructor selection table

        // For the dictionary constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            searchableSurface,
            dict,
            (
                const word& name,
                const objectRegistry& obj,
                const dictionary& dict
            ),
            (name, obj, dict)
        );

        //// For the Istream constructor
        //declareRunTimeSelectionTable
        //(
        //    autoPtr,
        //    searchableSurface,
        //    istream,
        //    (
        //        const objectRegistry& obj,
        //        Istream& is
        //    ),
        //    (obj, is)
        //);


        //- Class used for the read-construction of
        //  PtrLists of searchableSurface
        class iNew
        {
            const objectRegistry& obj_;

        public:

            iNew(const objectRegistry& obj)
            :
                obj_(obj)
            {}

            autoPtr<searchableSurface> operator()(Istream& is) const
            {
                word surfaceType(is);
                word name(is);
                dictionary dict(is);
                return searchableSurface::New(surfaceType, name, obj_, dict);
            }
        };


    // Constructors

        searchableSurface(const word& name);

        //- Clone
        virtual autoPtr<searchableSurface> clone() const
        {
            notImplemented("autoPtr<searchableSurface> clone() const");
            return autoPtr<searchableSurface>(NULL);
        }


    // Selectors

        //- Return a reference to the selected searchableSurface
        static autoPtr<searchableSurface> New
        (
            const word& surfaceType,
            const word& name,
            const objectRegistry& obj,
            const dictionary& dict
        );

        ////- Return a reference to the selected searchableSurface
        //static autoPtr<searchableSurface> New
        //(
        //    const word& surfaceType,
        //    const objectRegistry& obj,
        //    Istream& is
        //);


    // Destructor

        virtual ~searchableSurface();


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Calculate nearest point on surface. Returns
        //  - bool : any point found nearer than nearestDistSqr
        //  - label: relevant index in surface
        //  - point: actual nearest point found
        virtual pointIndexHit findNearest
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const = 0;

        //- Calculate nearest point on edge. Returns
        //  - bool : any point found nearer than nearestDistSqr
        //  - label: relevant index in surface
        //  - point: actual nearest point found
        virtual pointIndexHit findNearestOnEdge
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const = 0;

        //- Find nearest to segment. Returns
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
        ) const = 0;

        //- Find nearest intersection of line between start and end.
        virtual pointIndexHit findLine
        (
            const point& start,
            const point& end
        ) const = 0;

        //- Find any intersection of line between start and end.
        virtual pointIndexHit findLineAny
        (
            const point& start,
            const point& end
        ) const = 0;

        //- Determine type (inside/outside/mixed) for point. unknown if
        //  cannot be determined (e.g. non-manifold surface)
        virtual volumeType getVolumeType(const point&) const = 0;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
