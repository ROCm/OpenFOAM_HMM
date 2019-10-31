/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    A centralized collection of named coordinate systems.

Note
    Mixing normal constructors and the coordinateSystems::New constructor
    may yield unexpected results.

    \verbatim
    cat1
    {
        coordinateSystem
        {
            type   indirect;
            name   _10;
        }
        porosity    0.781;
        Darcy
        {
            d   d [0 -2 0 0 0]  (-1000 -1000 0.50753e+08);
            f   f [0 -1 0 0 0]  (-1000 -1000 12.83);
        }
    }
    \endverbatim

    For this to work correctly, the coordinateSystem constructor must be
    supplied with an objectRegistry as well as the dictionary.

SourceFiles
    coordinateSystems.C

\*---------------------------------------------------------------------------*/
#ifndef coordinateSystems_H
#define coordinateSystems_H

#include "regIOobject.H"
#include "PtrList.H"
#include "coordinateSystem.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coordinateSystems Declaration
\*---------------------------------------------------------------------------*/

class coordinateSystems
:
    public regIOobject,
    public PtrList<coordinateSystem>
{
    // Private Member Functions

        //- Read "coordinateSystems" or older "IOPtrList<coordinateSystem>"
        void readFromStream(const bool valid = true);

        //- Attempt read if MUST_READ.., or READ_IF_PRESENT and has header
        //  \return False if no read should have been attempted
        bool readObject(const IOobject& io);


        //- No copy construct
        coordinateSystems(const coordinateSystems&) = delete;

        //- No copy assignment
        void operator=(const coordinateSystems&) = delete;


public:

    //- Runtime type information
    TypeNameNoDebug("coordinateSystems");


    // Constructors

        //- Read construct from IOobject
        explicit coordinateSystems(const IOobject& io);

        //- Construct from IOobject and PtrList content
        coordinateSystems
        (
            const IOobject& io,
            const PtrList<coordinateSystem>& content
        );

        //- Construct from IOobject and transferring PtrList content
        coordinateSystems
        (
            const IOobject& io,
            PtrList<coordinateSystem>&& content
        );


    // Selectors

        //- Return previously registered or read construct from "constant"
        static const coordinateSystems& New(const objectRegistry& obr);


    // Member Functions

        //- Find and return indices for all matches
        labelList indices(const keyType& key) const;

        //- Find and return indices for all matches
        labelList indices(const wordRes& matcher) const;

        //- Find and return index for the first match, return -1 if not found
        label findIndex(const keyType& key) const;

        //- Find and return index for the first match, return -1 if not found
        label findIndex(const wordRes& matcher) const;

        //- Search if given key exists
        bool found(const keyType& key) const;

        //- Return reference to named coordinateSystem or FatalErrror
        const coordinateSystem& lookup(const word& name) const;

        //- Return pointer to named coordinateSystem or nullptr on error
        const coordinateSystem* lookupPtr(const word& name) const;

        //- A list of the coordinate-system names
        wordList names() const;

        //- A list of the coordinate-system names satisfying the input matcher
        wordList names(const keyType& key) const;

        //- A list of the coordinate-system names satisfying the input matcher
        wordList names(const wordRe& matcher) const;

        //- A list of the coordinate-system names satisfying the input matcher
        wordList names(const wordRes& matcher) const;

        //- Identical to names()
        inline wordList toc() const
        {
            return names();
        }


    // IO

        //- Write data
        bool writeData(Ostream& os) const;

        //- Write data
        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber ver,
            IOstream::compressionType,
            const bool valid = true
        ) const;


    // Housekeeping

        //- Identical to the indices() method (AUG-2018)
        labelList findIndices(const keyType& key) const
        {
            return this->indices(key);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
