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
    Foam::surfGroup

Description
    'Patch' on surface as subset of triSurface.

SourceFiles
    surfGroup.C

\*---------------------------------------------------------------------------*/

#ifndef surfGroup_H
#define surfGroup_H

#include "word.H"
#include "label.H"
#include "className.H"
#include "surfPatchIdentifier.H"
#include "autoPtr.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class surfGroup;

Istream& operator>>(Istream&, surfGroup&);
Ostream& operator<<(Ostream&, const surfGroup&);

/*---------------------------------------------------------------------------*\
                           Class surfGroup Declaration
\*---------------------------------------------------------------------------*/

class surfGroup
:
    public surfPatchIdentifier
{
    // Private data

        //- size of this group in the face list
        label size_;

        //- Start label of this group in the face list
        label start_;

public:

    //- Runtime type information
    ClassName("surfGroup");


    // Constructors

        //- Construct null
        surfGroup();

        //- Construct from components
        surfGroup
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const word& geometricType = word::null
        );

        //- Construct from Istream
        surfGroup(Istream& is, const label index);

        //- Construct from dictionary
        surfGroup
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Construct as copy
        surfGroup(const surfGroup&);

        //- Construct from another patch, resetting the index
        surfGroup
        (
            const surfGroup&,
            const label index
        );

        //- Return clone
        autoPtr<surfGroup> clone() const
        {
            notImplemented("autoPtr<surfGroup> clone() const");
            return autoPtr<surfGroup>(NULL);
        }

        static autoPtr<surfGroup> New(Istream& is)
        {
            word name(is);
            dictionary dict(is);

            return autoPtr<surfGroup>(new surfGroup(name, dict, 0));
        }


    // Member Functions

        //- Return start label of this patch in the face list
        label start() const
        {
            return start_;
        }

        //- Return start label of this patch in the face list
        label& start()
        {
            return start_;
        }

        //- Return size of this patch in the face list
        label size() const
        {
            return size_;
        }

        //- Return size of this patch in the face list
        label& size()
        {
            return size_;
        }

        //- Write
        void write(Ostream&) const;

        //- Write dictionary
        void writeDict(Ostream&) const;


    // Member Operators

        bool operator!=(const surfGroup&) const;

        //- compare.
        bool operator==(const surfGroup&) const;

    // IOstream Operators

        friend Istream& operator>>(Istream&, surfGroup&);
        friend Ostream& operator<<(Ostream&, const surfGroup&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
