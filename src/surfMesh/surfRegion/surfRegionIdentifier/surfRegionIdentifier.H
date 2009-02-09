/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::surfRegionIdentifier

Description
    An identifier for a region on a meshed surface.

    Similar in concept to a faceZone on the surface, but can also have a
    "geometricType" rather.  Despite the similarity to a 'patch' for
    volume meshes (with a "physicalType"), the region does not have any
    patch information per se.

SourceFiles
    surfRegionIdentifier.C

\*---------------------------------------------------------------------------*/

#ifndef surfRegionIdentifier_H
#define surfRegionIdentifier_H

#include "word.H"
#include "label.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

// Forward declaration of friend functions and operators

class surfRegionIdentifier;
Ostream& operator<<(Ostream&, const surfRegionIdentifier&);

/*---------------------------------------------------------------------------*\
                    Class surfRegionIdentifier Declaration
\*---------------------------------------------------------------------------*/

class surfRegionIdentifier
{
    // Private data

        //- Name of region
        word name_;

        //- Index of region in surface mesh
        label boundaryIndex_;

        //- Type name of region
        mutable word geometricType_;

public:

    // Constructors

        //- Construct null
        surfRegionIdentifier();

        //- Construct from components
        surfRegionIdentifier
        (
            const word& name,
            const label index,
            const word& geometricType = word::null

        );

        //- Construct from dictionary
        surfRegionIdentifier
        (
            const word& name,
            const dictionary&,
            const label index
        );

        //- Construct from another region identifier, resetting the index
        surfRegionIdentifier
        (
            const surfRegionIdentifier&,
            const label index
        );


    // Destructor

        virtual ~surfRegionIdentifier();


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return name for modification
        word& name()
        {
            return name_;
        }

        //- Return the geometric type of the region
        const word& geometricType() const
        {
            return geometricType_;
        }

        //- Return the geometric type of the region for modification
        word& geometricType()
        {
            return geometricType_;
        }

        //- Return the index of this region in the surface mesh
        label index() const
        {
            return boundaryIndex_;
        }

        //- Write surfRegionIdentifier as a dictionary
        void write(Ostream&) const;

        //- Write surfRegionIdentifier as a dictionary
//        void writeDict(Ostream&) const;


    // Member Operators

//        bool operator!=(const surfRegionIdentifier&) const;
//
//        //- compare.
//        bool operator==(const surfRegionIdentifier&) const;

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const surfRegionIdentifier&);
//        friend Istream& operator>>(Istream&, surfRegionIdentifier&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
