/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::geometricSurfacePatch

Description
    The geometricSurfacePatch is like patchIdentifier but for surfaces.
    Holds type, name and index.

SourceFiles
    geometricSurfacePatch.C

\*---------------------------------------------------------------------------*/

#ifndef geometricSurfacePatch_H
#define geometricSurfacePatch_H

#include "typeInfo.H"
#include "surfZoneIdentifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class dictionary;
class geometricSurfacePatch;

Istream& operator>>(Istream& is, geometricSurfacePatch& p);
Ostream& operator<<(Ostream& os, const geometricSurfacePatch& p);

/*---------------------------------------------------------------------------*\
                    Class geometricSurfacePatch Declaration
\*---------------------------------------------------------------------------*/

class geometricSurfacePatch
{
    // Private data

        //- Type name of patch
        word geometricType_;

        //- Name of patch
        word name_;

        //- Index of patch in boundary
        label index_;

public:

    // Public data

        //- The name for an 'empty' type
        static const word emptyType;


    //- Runtime type information
    ClassName("geometricSurfacePatch");


    // Constructors

        //- Construct null
        geometricSurfacePatch();

        //- Construct null with specified index
        explicit geometricSurfacePatch(const label index);

        //- Construct from components
        geometricSurfacePatch
        (
            const word& name,
            const label index,
            const word& geometricType = word::null
        );

        //- Construct from components
        geometricSurfacePatch
        (
            const word& geometricType,
            const word& name,
            const label index
        );

        //- Implicit conversion from surfZoneIdentifier, which has similar
        //- information but in a different order
        geometricSurfacePatch(const surfZoneIdentifier& ident);

        //- Construct from dictionary
        geometricSurfacePatch
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Construct from Istream
        geometricSurfacePatch(Istream& is, const label index);


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return name
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

        //- Return the index of this patch in the surface mesh
        label index() const
        {
            return index_;
        }

        //- Return the index of this patch in the surface mesh for modification
        label& index()
        {
            return index_;
        }

        //- Write
        void write(Ostream& os) const;

        //- Write dictionary
        void writeDict(Ostream& os) const;


    // Ostream Operator

        friend Istream& operator>>(Istream& is, geometricSurfacePatch& p);
        friend Ostream& operator<<(Ostream& os, const geometricSurfacePatch& p);
};


// Global Operators

//- Compare patches for equality
bool operator==(const geometricSurfacePatch& a, const geometricSurfacePatch& b);

//- Compare patches for inequality
bool operator!=(const geometricSurfacePatch& a, const geometricSurfacePatch& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
