/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::IOobjectList

Description
    List of IOobjects with searching and retrieving facilities.

SourceFiles
    IOobjectList.C

\*---------------------------------------------------------------------------*/

#ifndef IOobjectList_H
#define IOobjectList_H

#include "HashPtrTable.H"
#include "IOobject.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IOobjectList Declaration
\*---------------------------------------------------------------------------*/

class IOobjectList
:
    public HashPtrTable<IOobject>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const IOobjectList&) = delete;


public:

    // Constructors

        //- Construct given an initial estimate for the number of entries
        explicit IOobjectList(const label nIoObjects = 128);

        //- Construct from objectRegistry and instance path
        IOobjectList
        (
            const objectRegistry& db,
            const fileName& instance,
            const fileName& local = "",
            IOobject::readOption r = IOobject::MUST_READ,
            IOobject::writeOption w = IOobject::NO_WRITE,
            bool registerObject = true
        );

        //- Construct as copy
        IOobjectList(const IOobjectList& iolist);


    //- Destructor
    ~IOobjectList();


    // Member functions

        //- Add an IOobject to the list
        bool add(IOobject& io);

        //- Remove an IOobject from the list
        bool remove(IOobject& io);

        //- Lookup a given name and return IOobject ptr if found else nullptr
        IOobject* lookup(const word& name) const;

        //- The list of all IOobects with matching names
        IOobjectList lookup(const wordRe& matcher) const;

        //- The list of all IOobjects with matching names
        IOobjectList lookup(const wordReList& matcher) const;

        //- The list of all IOobjects with the given class name
        IOobjectList lookupClass(const word& clsName) const;


        //- A list of names of the IOobjects
        wordList names() const;

        //- The names of IOobjects with the given class name
        wordList names(const word& clsName) const;

        //- The names of IOobjects with the given class name that also
        //  have a name satisfying the input matcher
        wordList names(const word& clsName, const wordRe& matcher) const;

        //- The names of IOobjects with the given class name that also
        //  have a name satisfying the input matcher
        wordList names(const word& clsName, const wordReList& matcher) const;


        //- A sorted list of names of the IOobjects
        wordList sortedNames() const;

        //- The sorted names of IOobjects with the given class name
        wordList sortedNames(const word& clsName) const;

        //- The sorted names of IOobjects with the given class name that also
        //  have a name satisfying the input matcher
        wordList sortedNames(const word& clsName, const wordRe& matcher) const;

        //- The sorted names of IOobjects with the given class name that also
        //  have a name satisfying the input matcher
        wordList sortedNames
        (
            const word& clsName,
            const wordReList& matcher
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
