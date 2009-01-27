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
    Foam::surfRegion

Description
    A region on a meshed surface.
    Similar in concept to a faceZone on the surface.

SourceFiles
    surfRegion.C

\*---------------------------------------------------------------------------*/

#ifndef surfRegion_H
#define surfRegion_H

#include "word.H"
#include "label.H"
#include "className.H"
#include "surfRegionIdentifier.H"
#include "autoPtr.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class surfRegion;

Istream& operator>>(Istream&, surfRegion&);
Ostream& operator<<(Ostream&, const surfRegion&);

/*---------------------------------------------------------------------------*\
                           Class surfRegion Declaration
\*---------------------------------------------------------------------------*/

class surfRegion
:
    public surfRegionIdentifier
{
    // Private data

        //- size of this group in the face list
        label size_;

        //- Start label of this group in the face list
        label start_;

public:

    //- Runtime type information
    ClassName("surfRegion");


    // Constructors

        //- Construct null
        surfRegion();

        //- Construct from components
        surfRegion
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const word& geometricType = word::null
        );

        //- Construct from Istream
        surfRegion(Istream& is, const label index);

        //- Construct from dictionary
        surfRegion
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Construct as copy
        surfRegion(const surfRegion&);

        //- Construct from another region, resetting the index
        surfRegion(const surfRegion&, const label index);

        //- Return clone
        autoPtr<surfRegion> clone() const
        {
            notImplemented("autoPtr<surfRegion> clone() const");
            return autoPtr<surfRegion>(NULL);
        }

        static autoPtr<surfRegion> New(Istream& is)
        {
            word name(is);
            dictionary dict(is);

            return autoPtr<surfRegion>(new surfRegion(name, dict, 0));
        }


    // Member Functions

        //- Return start label of this region in the face list
        label start() const
        {
            return start_;
        }

        //- Return start label of this region in the face list
        label& start()
        {
            return start_;
        }

        //- Return size of this region in the face list
        label size() const
        {
            return size_;
        }

        //- Return size of this region in the face list
        label& size()
        {
            return size_;
        }

        //- Write
        void write(Ostream&) const;

        //- Write dictionary
        void writeDict(Ostream&) const;


    // Member Operators

        bool operator!=(const surfRegion&) const;

        //- compare.
        bool operator==(const surfRegion&) const;

    // IOstream Operators

        friend Istream& operator>>(Istream&, surfRegion&);
        friend Ostream& operator<<(Ostream&, const surfRegion&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
