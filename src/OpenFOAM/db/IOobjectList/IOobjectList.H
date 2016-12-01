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
        IOobjectList(const IOobjectList&);


    //- Destructor
    ~IOobjectList();


    // Member functions

        //- Add an IOobject to the list
        bool add(IOobject&);

        //- Remove an IOobject from the list
        bool remove(IOobject&);

        //- Lookup a given name and return IOobject ptr if found else nullptr
        IOobject* lookup(const word& name) const;

        //- Return the list for all IOobects whose name matches name
        IOobjectList lookup(const wordRe& name) const;

        //- Return the list for all IOobects whose name matches name
        IOobjectList lookup(const wordReList& patterns) const;

        //- Return the list for all IOobjects of a given class
        IOobjectList lookupClass(const word& className) const;


        //- A list of names of the IOobjects
        wordList names() const;

        //- A list of names of IOobjects of the given class
        wordList names(const word& className) const;

        //- A list of names of IOobjects of the given class,
        //  and that also satisfy the input matcher
        wordList names(const word& className, const wordRe&) const;

        //- A list of names of IOobjects of the given class,
        //  and that also satisfy the input matchers
        wordList names(const word& className, const wordReList&) const;


        //- A sorted list of names of the IOobjects
        wordList sortedNames() const;

        //- A sorted list of names of IOobjects of given class
        wordList sortedNames(const word& className) const;

        //- A sorted list of names of IOobjects of the given class,
        //  and that also satisfy the input matcher
        wordList sortedNames(const word& className, const wordRe&) const;

        //- A sorted list of names of IOobjects of the given class,
        //  and that also satisfy the input matchers
        wordList sortedNames(const word& className, const wordReList&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
