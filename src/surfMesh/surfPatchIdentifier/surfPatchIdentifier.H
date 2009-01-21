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
    Foam::surfPatchIdentifier

Description
    Like patchIdentifier but for surfaces with "geometricType" rather
    than "physicalType".

SourceFiles
    surfPatchIdentifier.C

\*---------------------------------------------------------------------------*/

#ifndef surfPatchIdentifier_H
#define surfPatchIdentifier_H

#include "word.H"
#include "label.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

// Forward declaration of friend functions and operators

class surfPatchIdentifier;
Ostream& operator<<(Ostream&, const surfPatchIdentifier&);

/*---------------------------------------------------------------------------*\
                    Class surfPatchIdentifier Declaration
\*---------------------------------------------------------------------------*/

class surfPatchIdentifier
{
    // Private data

        //- Name of patch
        word name_;

        //- Index of patch in boundary
        label boundaryIndex_;

        //- Type name of patch
        mutable word geometricType_;

public:

    // Constructors

        //- Construct null
        surfPatchIdentifier();

        //- Construct from components
        surfPatchIdentifier
        (
            const word& name,
            const label index,
            const word& geometricType = word::null

        );

        //- Construct from dictionary
        surfPatchIdentifier
        (
            const word& name,
            const dictionary&,
            const label index
        );

        //- Construct from another patch, resetting the index
        surfPatchIdentifier
        (
            const surfPatchIdentifier&,
            const label index
        );


    // Destructor

        virtual ~surfPatchIdentifier();


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

        //- Return the geometric type of the patch
        const word& geometricType() const
        {
            return geometricType_;
        }

        //- Return the geometric type of the patch for modification
        word& geometricType()
        {
            return geometricType_;
        }

        //- Return the index of this patch in the boundaryMesh
        label index() const
        {
            return boundaryIndex_;
        }

        //- Write surfPatchIdentifier as a dictionary
        void write(Ostream&) const;

        //- Write surfPatchIdentifier as a dictionary
//        void writeDict(Ostream&) const;


    // Member Operators

//        bool operator!=(const surfPatchIdentifier&) const;
//
//        //- compare.
//        bool operator==(const surfPatchIdentifier&) const;

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const surfPatchIdentifier&);
//        friend Istream& operator>>(Istream&, surfPatchIdentifier&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
