/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
    IOobjectListTemplates.C

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

/*---------------------------------------------------------------------------*\
                        Class IOobjectList Declaration
\*---------------------------------------------------------------------------*/

class IOobjectList
:
    public HashPtrTable<IOobject>
{
    // Private Member Functions

        //- Check name consistency on all processors
        //
        //  With syncPar = true, check that object names are the same on
        //  all processors. Trigger FatalError if not.
        //
        //  When syncPar is used, the object names are sorted as a side-effect,
        //  since this is required for consistent ordering across processors.
        static bool checkNames(wordList& masterNames, const bool syncPar);

        //- Combine names from all processors and sort
        static void syncNames(wordList& objNames);

        //- Templated implementation for classes()
        template<class MatchPredicate>
        static HashTable<wordHashSet> classesImpl
        (
            const IOobjectList& list,
            const MatchPredicate& matchName
        );

        //- Templated implementation for count()
        //  The number of items with a matching class
        template<class MatchPredicate1, class MatchPredicate2>
        static label countImpl
        (
            const IOobjectList& list,
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        );

        //- Templated implementation for count()
        template<class Type, class MatchPredicate>
        static label countTypeImpl
        (
            const IOobjectList& list,
            const MatchPredicate& matchName
        );

        //- Templated implementation for names(), sortedNames()
        template<class MatchPredicate1, class MatchPredicate2>
        static wordList namesImpl
        (
            const IOobjectList& list,
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName,
            const bool doSort
        );

        //- Templated implementation for names(), sortedNames()
        template<class Type, class MatchPredicate>
        static wordList namesTypeImpl
        (
            const IOobjectList& list,
            const MatchPredicate& matchName,
            const bool doSort
        );

        //- Templated implementation for lookup()
        template<class MatchPredicate>
        static IOobjectList lookupImpl
        (
            const IOobjectList& list,
            const MatchPredicate& matchName
        );

        //- Templated implementation for lookupClass()
        template<class MatchPredicate1, class MatchPredicate2>
        static IOobjectList lookupClassImpl
        (
            const IOobjectList& list,
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        );

        //- Templated implementation for lookupClass()
        template<class Type, class MatchPredicate>
        static IOobjectList lookupClassTypeImpl
        (
            const IOobjectList& list,
            const MatchPredicate& matchName
        );


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


    // Member Functions

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


    // Lookup single item

        //- Return const pointer to the object found by name
        //  \return IOobject ptr if found else nullptr
        const IOobject* cfindObject(const word& objName) const;

        //- Return const pointer to the object found by name that also
        //- has headerClassName == Type::typeName
        //  \return IOobject ptr if found and the correct type, else nullptr
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type>
        const IOobject* cfindObject(const word& objName) const;

        //- Return const pointer to the object found by name
        //  \return IOobject ptr if found else nullptr
        const IOobject* findObject(const word& objName) const;

        //- Return const pointer to the object found by name that also
        //- has headerClassName == Type::typeName
        //  \return IOobject ptr if found and the correct type, else nullptr
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type>
        const IOobject* findObject(const word& objName) const;

        //- Return non-const pointer to the object found by name
        //  \return IOobject ptr if found else nullptr
        IOobject* findObject(const word& objName);

        //- Return non-const pointer to the object found by name that also
        //- has headerClassName == Type::typeName
        //  \return IOobject ptr if found and the correct type, else nullptr
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type>
        IOobject* findObject(const word& objName);

        //- Return non-const pointer to the object found by name,
        //- using a const-cast to have it behave like a mutable.
        //  Exercise caution when using.
        //  \return IOobject ptr if found else nullptr
        IOobject* getObject(const word& objName) const;

        //- Return non-const pointer to the object found by name that also
        //- has headerClassName == Type::typeName,
        //- using a const-cast to have it behave like a mutable.
        //  Exercise caution when using.
        //  \return IOobject ptr if found and the correct type, else nullptr
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type>
        IOobject* getObject(const word& objName) const;


    // Lookup multiple items

        //- The list of IOobjects that have a matching object name.
        template<class MatchPredicate>
        IOobjectList lookup(const MatchPredicate& matchName) const;

        //- The list of IOobjects with the given headerClassName
        IOobjectList lookupClass(const char* clsName) const;

        //- The list of IOobjects with matching headerClassName
        template<class MatchPredicate>
        IOobjectList lookupClass(const MatchPredicate& matchClass) const;

        //- The list of IOobjects with matching headerClassName
        //- that also have a matching object name.
        template<class MatchPredicate1, class MatchPredicate2>
        IOobjectList lookupClass
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The list of IOobjects with headerClassName == Type::typeName
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type>
        IOobjectList lookupClass() const;

        //- The list of IOobjects with headerClassName == Type::typeName
        //- that also have a matching object name.
        template<class Type, class MatchPredicate>
        IOobjectList lookupClass(const MatchPredicate& matchName) const;


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
        //
        //  \return HashTable of class-names for keys and wordHashSet of
        //      of object-names for the values.
        HashTable<wordHashSet> classes() const;

        //- A summary hash of classes used and their associated object names,
        //- restricted to objects that have a matching object name.
        template<class MatchPredicate>
        HashTable<wordHashSet> classes(const MatchPredicate& matchName) const;


    // Number of items

        //- The number of objects of the given headerClassName
        //  \note uses the class type() method
        label count(const char* clsName) const;

        //- The number of objects of the given headerClassName
        template<class MatchPredicate>
        label count
        (
            const MatchPredicate& matchClass
        ) const;

        //- The number of objects of the given headerClassName
        //- that also have a matching object name.
        template<class MatchPredicate1, class MatchPredicate2>
        label count
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The number of objects with headerClassName == Type::typeName
        template<class Type>
        label count() const;

        //- The number of objects with headerClassName == Type::typeName
        //- that also have a matching object name.
        //
        //  \note If \a Type is \c void, no headerClassName check is used
        //      (always true).
        template<class Type, class MatchPredicate>
        label count(const MatchPredicate& matchName) const;


    // Summary of names

        //- The names of the IOobjects
        wordList names() const;

        //- The names of the IOobjects.
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        wordList names(const bool syncPar) const;

        //- The names of IOobjects with the given headerClassName
        wordList names(const char* clsName) const;

        //- The names of the IOobjects with the given headerClassName
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        wordList names(const char* clsName, const bool syncPar) const;

        //- The names of IOobjects with the given headerClassName
        template<class MatchPredicate>
        wordList names(const MatchPredicate& matchClass) const;

        //- The names of the IOobjects with the given headerClassName
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        template<class MatchPredicate>
        wordList names
        (
            const MatchPredicate& matchClass,
            const bool syncPar
        ) const;

        //- The names of IOobjects with the given headerClassName
        //- that also have a matching object name.
        template<class MatchPredicate1, class MatchPredicate2>
        wordList names
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The names of the IOobjects with the given headerClassName
        //- that also have a matching object name.
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        template<class MatchPredicate1, class MatchPredicate2>
        wordList names
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName,
            const bool syncPar
        ) const;


        //- The names of objects with headerClassName == Type::typeName
        template<class Type>
        wordList names() const;

        //- The names of objects with headerClassName == Type::typeName
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        template<class Type>
        wordList names(bool syncPar) const;

        //- The names of objects with headerClassName == Type::typeName
        //- that also have a matching object name.
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        template<class Type, class MatchPredicate>
        wordList names(const MatchPredicate& matchName) const;

        //- The names of objects with headerClassName == Type::typeName
        //- that also have a matching object name.
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        template<class Type, class MatchPredicate>
        wordList names
        (
            const MatchPredicate& matchName,
            const bool syncPar
        ) const;


    // Summary of names (sorted)

        //- The sorted names of the IOobjects
        wordList sortedNames() const;

        //- The sorted names of the IOobjects.
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        wordList sortedNames(const bool syncPar) const;

        //- The sorted names of IOobjects with the given headerClassName
        wordList sortedNames(const char* clsName) const;

        //- The sorted names of the IOobjects with the given headerClassName
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        wordList sortedNames(const char* clsName, const bool syncPar) const;

        //- The sorted names of IOobjects with the given headerClassName
        template<class MatchPredicate>
        wordList sortedNames(const MatchPredicate& matchClass) const;

        //- The sorted names of the IOobjects with the given headerClassName
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        template<class MatchPredicate>
        wordList sortedNames
        (
            const MatchPredicate& matchClass,
            const bool syncPar
        ) const;

        //- The sorted names of IOobjects with the given headerClassName
        //- that also have a matching object name.
        template<class MatchPredicate1, class MatchPredicate2>
        wordList sortedNames
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The sorted names of the IOobjects with the given headerClassName
        //- that also have a matching object name.
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        template<class MatchPredicate1, class MatchPredicate2>
        wordList sortedNames
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName,
            const bool syncPar
        ) const;


        //- The sorted names of objects with headerClassName == Type::typeName
        template<class Type>
        wordList sortedNames() const;

        //- The sorted names of objects with headerClassName == Type::typeName
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        template<class Type>
        wordList sortedNames(bool syncPar) const;

        //- The sorted names of objects with headerClassName == Type::typeName
        //- that also have a matching object name.
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        template<class Type, class MatchPredicate>
        wordList sortedNames(const MatchPredicate& matchName) const;

        //- The sorted names of objects with headerClassName == Type::typeName
        //- that also have a matching object name.
        //  With syncPar = true, sorts the names and triggers FatalError
        //  if the names are not consistent on all processors.
        template<class Type, class MatchPredicate>
        wordList sortedNames
        (
            const MatchPredicate& matchName,
            const bool syncPar
        ) const;


    // Edit

        //- Filter to retain or prune given classes
        //  \return The number of items changed (removed)
        template<class UnaryPredicate>
        label filterClasses
        (
            const UnaryPredicate& pred,
            const bool pruning = false
        );

        //- Filter to retain or prune given object names
        //  \return The number of items changed (removed)
        template<class UnaryPredicate>
        label filterObjects
        (
            const UnaryPredicate& pred,
            const bool pruning = false
        );

        //- Remove objects with names ending with "_0" (restart fields)
        //  \return The number of items changed (removed)
        label prune_0();


    // Parallel

        //- The sorted names of all objects (synchronised across processors)
        wordList allNames() const;

        //- The sorted names of all objects (synchronised across processors)
        //- with headerClassName == Type::typeName
        template<class Type>
        wordList allNames() const;

        //- Verify that object names are synchronised across processors
        //  Triggers FatalError if the names are not consistent on all
        //  processors.
        bool checkNames(const bool syncPar = true) const;


    // Member Operators

        //- No copy assignment
        void operator=(const IOobjectList&) = delete;

        //- Move assignment
        void operator=(IOobjectList&& list);


    // Housekeeping

        //- Deprecated(2018-11) Locate an object by name (c-string).
        //- Disambiguated from multiple-lookup version by calling parameter.
        //  \deprecated(2018-11) use findObject() for non-ambiguous resolution
        IOobject* lookup(const char* objName) const
        {
            return getObject(objName);
        }

        //- Deprecated(2018-11) Locate an object by name (const word&).
        //- Disambiguated from multiple-lookup version by calling parameter.
        //  \deprecated(2018-11) use findObject() for non-ambiguous resolution
        IOobject* lookup(const word& objName) const
        {
            return getObject(objName);
        }
};


// Ostream operator
Ostream& operator<<(Ostream& os, const IOobjectList& list);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOobjectListTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
