/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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
#include "HashSet.H"
#include "IOobject.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class IOobjectList;
Ostream& operator<<(Ostream& os, const IOobjectList& list);


/*---------------------------------------------------------------------------*\
                        Class IOobjectList Declaration
\*---------------------------------------------------------------------------*/

class IOobjectList
:
    public HashPtrTable<IOobject>
{
public:

    // Constructors

        //- Construct null with default (128) table capacity
        IOobjectList();

        //- Construct given initial table capacity
        explicit IOobjectList(const label nObjects);

        //- Copy construct
        IOobjectList(const IOobjectList& list);

        //- Move construct
        IOobjectList(IOobjectList&& list);

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


    //- Destructor
    ~IOobjectList() = default;


    // Member functions

    // Basic methods

        //- Add IOobject to the list
        bool add(autoPtr<IOobject>& objectPtr);

        //- Add IOobject to the list
        bool add(autoPtr<IOobject>&& objectPtr);

        //- Copy append objects from other to this list, but do not overwrite
        //- existing keys.
        //
        //  \return number of items added
        label append(const IOobjectList& other);

        //- Move append objects from other to this list, but do not overwrite
        //- existing keys.
        //  After calling this, the other parameter will contains any items
        //  that could not be moved.
        //
        //  \return number of items added
        label append(IOobjectList&& other);

        //- Remove IOobject from the list, by name or by iterator.
        //
        //  \return autoPtr<IOobject>
        using HashPtrTable<IOobject>::remove;

        //- Remove IOobject from the list.
        //
        //  \return True if object was removed
        bool remove(const IOobject& io);


    // Lookup

        //- Lookup a given name and return IOobject ptr if found else nullptr
        IOobject* lookup(const word& name) const;

        //- The list of all IOobjects with matching names
        IOobjectList lookup(const wordRe& matcher) const;

        //- The list of all IOobjects with matching names
        IOobjectList lookup(const wordRes& matcher) const;

        //- The list of all IOobjects with matching names
        IOobjectList lookup(const wordHashSet& matcher) const;

        //- The list of all IOobjects with the given class name
        IOobjectList lookupClass(const word& clsName) const;


    // Summary of classes

        //- A summary hash of classes used and their associated object names.
        //  The HashTable representation allows us to leverage various
        //  HashTable methods.
        //  This hashed summary view can be useful when querying particular
        //  aspects. For example,
        //
        //  \code
        //  IOobjectList objects(runTime, runTime.timeName());
        //  HashTable<wordHashSet> classes = objects.classes();
        //
        //  // How many volScalarField?
        //  word checkType = volScalarField::typeName;
        //
        //  Info<< checkType << "="
        //      << (classes.found(checkType) ? classes[checkType].size() : 0)
        //      << nl;
        //  \endcode
        //  Using the two-parameter HashTable::lookup method lets us avoid
        //  the \c '?' ternary, but still looks fairly ugly:
        //  \code
        //  Info<< checkType << "="
        //      << classes.lookup(checkType, wordHashSet()).size() << nl;
        //  \endcode
        //
        //  If we have non-const access to the hash table, and don't mind
        //  incidentally creating empty entries,
        //  we can use the HashTable::operator() directly:
        //  \code
        //  Info<< checkType << "=" << classes(checkType).size() << nl;
        //  \endcode
        //
        //  Of course, for a single query it would have been easier
        //  and simpler to have used a direct query of the names:
        //  \code
        //  Info<< checkType << "=" << objects.names(checkType).size() << nl;
        //  \endcode
        //
        //  The summary hash, however, becomes most useful when reducing
        //  the objects in consideration to a particular subset. For example,
        //  \code
        //  const wordHashSet interestingTypes
        //  {
        //      volScalarField::typeName,
        //      volVectorField::typeName
        //  };
        //
        //  classes.retain(interestingTypes);
        //  \endcode
        //  Or do just the opposite:
        //  \code
        //  classes.erase(unsupportedTypes);
        //  \endcode
        //  This also works with a hashedWordList, since it provides the
        //  expected '()' operator. But in this case the more general
        //  HashTable::filterKeys is required:
        //  \code
        //  const hashedWordList interestingTypes
        //  {
        //      volScalarField::typeName,
        //      volVectorField::typeName
        //  };
        //
        //  classes.filterKeys(interestingTypes);
        //  \endcode
        //
        //  Of course, there are many other ways to use and manipulate the
        //  summary information.
        HashTable<wordHashSet> classes() const;

        //- A summary hash of classes used and their associated object names
        //  restricted to objects with names that satisfy the input matcher
        HashTable<wordHashSet> classes(const wordRe& matcher) const;

        //- A summary hash of classes used and their associated object names
        //  restricted to objects with names that satisfy the input matcher
        HashTable<wordHashSet> classes(const wordRes& matcher) const;

        //- A summary hash of classes used and their associated object names
        //  restricted to objects with names that satisfy the input matcher
        HashTable<wordHashSet> classes(const wordHashSet& matcher) const;


    // Summary of names

        //- A list of names of the IOobjects
        wordList names() const;

        //- The names of IOobjects with the given class name
        wordList names(const word& clsName) const;

        //- The names of IOobjects with the given class name that also
        //- have a name satisfying the input matcher
        wordList names(const word& clsName, const wordRe& matcher) const;

        //- The names of IOobjects with the given class name that also
        //- have a name satisfying the input matcher
        wordList names(const word& clsName, const wordRes& matcher) const;

        //- The names of IOobjects with the given class name that also
        //- have a name satisfying the input matcher
        wordList names(const word& clsName, const wordHashSet& matcher) const;


        //- A sorted list of names of the IOobjects.
        //  With syncPar = true, triggers FatalError if the names are
        //  not consistent on all processors.
        wordList names(const bool syncPar) const;

        //- A sorted list of names of the IOobjects
        //  With syncPar = true, triggers FatalError if the names are
        //  not consistent on all processors.
        wordList names(const word& clsName, const bool syncPar) const;

        //- The sorted names of IOobjects with the given class name that also
        //- have a name satisfying the input matcher
        //  With syncPar = true, triggers FatalError if the names are
        //  not consistent on all processors.
        wordList names
        (
            const word& clsName,
            const wordRe& matcher,
            const bool syncPar
        ) const;

        //- The sorted names of IOobjects with the given class name that also
        //- have a name satisfying the input matcher
        //  With syncPar = true, triggers FatalError if the names are
        //  not consistent on all processors.
        wordList names
        (
            const word& clsName,
            const wordRes& matcher,
            const bool syncPar
        ) const;

        //- The sorted names of IOobjects with the given class name that also
        //- have a name satisfying the input matcher
        //  With syncPar = true, triggers FatalError if the names are
        //  not consistent on all processors.
        wordList names
        (
            const word& cls,
            const wordHashSet& matcher,
            const bool syncPar
        ) const;


    // Summary of names (sorted)

        //- A sorted list of names of the IOobjects
        wordList sortedNames() const;

        //- The sorted names of IOobjects with the given class name
        wordList sortedNames(const word& clsName) const;

        //- The sorted names of IOobjects with the given class name that also
        //- have a name satisfying the input matcher
        wordList sortedNames(const word& clsName, const wordRe& matcher) const;

        //- The sorted names of IOobjects with the given class name that also
        //- have a name satisfying the input matcher
        wordList sortedNames(const word& clsName, const wordRes& matcher) const;

        //- The sorted names of IOobjects with the given class name that also
        //- have a name satisfying the input matcher
        wordList sortedNames
        (
            const word& clsName,
            const wordHashSet& matcher
        ) const;


    // Member Operators

        //- No copy assignment
        void operator=(const IOobjectList&) = delete;

        //- Move assignment
        void operator=(IOobjectList&& list);


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const IOobjectList& list);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
