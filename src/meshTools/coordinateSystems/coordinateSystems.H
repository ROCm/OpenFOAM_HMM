/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::coordinateSystems

Description
    Provides a centralized coordinateSystem collection.

Note
    Mixing normal constructors and the coordinateSystems::New constructor
    may yield unexpected results.

    \verbatim
        1
        (
        cat1
        {
            coordinateSystem  system_10;
            porosity        0.781;
            Darcy
            {
                d   d [0 -2 0 0 0]  (-1000 -1000 0.50753e+08);
                f   f [0 -1 0 0 0]  (-1000 -1000 12.83);
            }
        }
        )
    \endverbatim

    For this to work correctly, the coordinateSystem constructor must be
    supplied with both a dictionary and an objectRegistry.

SourceFiles
    coordinateSystems.C

\*---------------------------------------------------------------------------*/
#ifndef coordinateSystems_H
#define coordinateSystems_H

#include "coordinateSystem.H"
#include "IOPtrList.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coordinateSystems Declaration
\*---------------------------------------------------------------------------*/

class coordinateSystems
:
    public IOPtrList<coordinateSystem>
{
    // Private Member Functions

        //- No copy construct
        coordinateSystems(const coordinateSystems&) = delete;

        //- No copy assignment
        void operator=(const coordinateSystems&) = delete;


public:

    //- Runtime type information
    TypeName("coordinateSystems");

    // Constructors

        //- Read construct from IOobject
        explicit coordinateSystems(const IOobject& io);

        //- Construct from IOobject and a PtrList
        coordinateSystems
        (
            const IOobject& io,
            const PtrList<coordinateSystem>& content
        );

        //- Construct from IOobject and transferring the PtrList content
        coordinateSystems
        (
            const IOobject& io,
            PtrList<coordinateSystem>&& content
        );


    // Selectors

        //- Return previously registered or read construct from "constant"
        static const coordinateSystems& New(const objectRegistry&);


    // Member Functions

        //- Find and return indices for all matches
        labelList findIndices(const keyType& key) const;

        //- Find and return index for the first match, return -1 if not found
        label findIndex(const keyType& key) const;

        //- Search for given key
        bool found(const keyType& key) const;

        //- A list of the coordinate-system names
        wordList names() const;

        //- A list of the coordinate-system names satisfying the input matcher
        wordList names(const keyType& matcher) const;

        //- A list of the coordinate-system names satisfying the input matcher
        wordList names(const wordRe& matcher) const;

        //- A list of the coordinate-system names satisfying the input matcher
        wordList names(const wordRes& matcher) const;

        //- Identical to names()
        inline wordList toc() const
        {
            return names();
        }

        //- Write data
        bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
