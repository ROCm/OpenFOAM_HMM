/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::surfZoneIdentifier

Description
    Identifies a surface patch/zone by name, patch index and geometricType.

SeeAlso
    patchIdentifier

SourceFiles
    surfZoneIdentifier.C

\*---------------------------------------------------------------------------*/

#ifndef surfZoneIdentifier_H
#define surfZoneIdentifier_H

#include "word.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                     Class surfZoneIdentifier Declaration
\*---------------------------------------------------------------------------*/

class surfZoneIdentifier
{
    // Private data

        //- Name of zone
        word name_;

        //- Index of zone in surface mesh
        label index_;

        //- Type name of zone
        mutable word geometricType_;


public:

    // Public data

        //- The name for an 'empty' type
        static const word emptyType;


    // Generated Methods

        //- Copy construct
        surfZoneIdentifier(const surfZoneIdentifier&) = default;

        //- Copy assignment
        surfZoneIdentifier& operator=(const surfZoneIdentifier&) = default;


    // Constructors

        //- Construct null, with index zero
        surfZoneIdentifier();

        //- Construct null with specified index
        explicit surfZoneIdentifier(const label index);

        //- Construct from components
        surfZoneIdentifier
        (
            const word& name,
            const label index,
            const word& geometricType = word::null
        );

        //- Construct from dictionary
        surfZoneIdentifier
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Copy construct from another zone identifier, resetting the index
        surfZoneIdentifier
        (
            const surfZoneIdentifier& p,
            const label index
        );


    //- Destructor
    virtual ~surfZoneIdentifier() = default;


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

        //- Return the geometric type of the patch/zone
        const word& geometricType() const
        {
            return geometricType_;
        }

        //- Return the geometric type of the patch/zone for modification
        word& geometricType()
        {
            return geometricType_;
        }

        //- Return the index of this patch/zone in the surface mesh
        label index() const
        {
            return index_;
        }

        //- Return the index of this patch/zone for modification
        label& index()
        {
            return index_;
        }


        //- Write identifier as a dictionary
        void write(Ostream& os) const;
};


// Global Operators

// Ostream Operator

//- Read name/type
Istream& operator>>(Istream& is, surfZoneIdentifier& obj);

//- Write name/type
Ostream& operator<<(Ostream& os, const surfZoneIdentifier& obj);

//- Compare zone indentifiers for equality
bool operator==(const surfZoneIdentifier& a, const surfZoneIdentifier& b);

//- Compare zone indentifiers for inequality
bool operator!=(const surfZoneIdentifier& a, const surfZoneIdentifier& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
