/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::objectRegistry

Description
    Registry of regIOobjects

SourceFiles
    objectRegistry.C
    objectRegistryTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef objectRegistry_H
#define objectRegistry_H

#include "HashTable.H"
#include "HashSet.H"
#include "regIOobject.H"
#include "wordRes.H"

// Historically included by objectRegistryTemplates (until NOV-2018),
// but not used by objectRegistry directly.
// Leave here for now to avoid a missing include in other bits of code.
#include "stringListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class objectRegistry Declaration
\*---------------------------------------------------------------------------*/

class objectRegistry
:
    public regIOobject,
    public HashTable<regIOobject*>
{
    // Private Data

        //- Master time objectRegistry
        const Time& time_;

        //- Parent objectRegistry
        const objectRegistry& parent_;

        //- Local directory path of this objectRegistry relative to time
        fileName dbDir_;

        //- Current event
        mutable label event_;


    // Private Member Functions

        //- Is the objectRegistry parent_ different from time_
        //  Used to terminate searching within the ancestors
        bool parentNotTime() const noexcept;

        //- Templated implementation for count()
        //  The number of items with a matching class
        template<class MatchPredicate1, class MatchPredicate2>
        static label countImpl
        (
            const objectRegistry& list,
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        );

        //- Templated implementation for count()
        //  The number of items with a matching class
        template<class Type, class MatchPredicate>
        static label countTypeImpl
        (
            const objectRegistry& list,
            const MatchPredicate& matchName
        );

        //- Templated implementation for classes()
        template<class MatchPredicate>
        static HashTable<wordHashSet> classesImpl
        (
            const objectRegistry& list,
            const MatchPredicate& matchName
        );

        //- Templated implementation for names(), sortedNames()
        template<class MatchPredicate1, class MatchPredicate2>
        static wordList namesImpl
        (
            const objectRegistry& list,
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName,
            const bool doSort
        );

        //- Templated implementation for names(), sortedNames()
        template<class Type, class MatchPredicate>
        static wordList namesTypeImpl
        (
            const objectRegistry& list,
            const MatchPredicate& matchName,
            const bool doSort
        );


        //- No copy construct
        objectRegistry(const objectRegistry&) = delete;

        //- No copy assignment
        void operator=(const objectRegistry&) = delete;


public:

    //- Declare type name for this IOobject
    TypeName("objectRegistry");


    // Constructors

        //- Construct the time objectRegistry,
        //- with estimated table capacity (default: 128)
        explicit objectRegistry(const Time& db, const label nObjects=128);

        //- Construct sub-registry given an IObject to describe the registry,
        //- with estimated table capacity (default: 128)
        explicit objectRegistry(const IOobject& io, const label nObjects=128);


    //- Destructor, with checkOut() for all objects that are ownedByRegistry
    virtual ~objectRegistry();


    // Member Functions

    // Access

        //- Return the object registry
        const objectRegistry& thisDb() const noexcept
        {
            return *this;
        }

        //- Return the parent objectRegistry
        const objectRegistry& parent() const noexcept
        {
            return parent_;
        }

        //- Return time registry
        const Time& time() const noexcept
        {
            return time_;
        }

        //- True if the registry is Time
        bool isTimeDb() const noexcept;

        //- Local directory path of this objectRegistry relative to the time
        virtual const fileName& dbDir() const
        {
            return dbDir_;
        }


    // Summary of classes

        //- A summary hash of classes used and their associated object names.
        //  Behaviour and usage as per IOobjectList::classes
        HashTable<wordHashSet> classes() const;

        //- A summary hash of classes used and their associated object names,
        //- restricted to objects that have a matching object name.
        template<class MatchPredicate>
        HashTable<wordHashSet> classes(const MatchPredicate& matchName) const;


    // Number of items

        //- The number of objects of the given class name
        //  \note uses the class type() method
        label count(const char* clsName) const;

        //- The number of objects of the given class name
        //  \note uses the class type() method
        template<class MatchPredicate>
        label count(const MatchPredicate& matchClass) const;

        //- The number of objects of the given class name
        //  \note uses the class type() method
        template<class MatchPredicate1, class MatchPredicate2>
        label count
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The names of objects with a class satisfying \c isA\<Type\>
        //
        //  \param strict use \c isType\<Type\> instead of \c isA\<Type\>
        //
        //  \note The values of \c count\<Type\>() and \c count(Type::typeName)
        //      may be inconsistent, since they use different mechanisms for
        //      testing the class type.
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type>
        label count(const bool strict = false) const;

        //- The names of objects with a class satisfying \c isA\<Type\>
        //- that also have a matching object name.
        //
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type, class MatchPredicate>
        label count(const MatchPredicate& matchName) const;


    // Summary of names

        //- The names of all objects
        wordList names() const;

        //- The names of objects with the given class name.
        //  \note uses the class type() method
        wordList names(const char* clsName) const;

        //- The names of objects with a matching class name
        //  \note uses the class type() method
        template<class MatchPredicate>
        wordList names(const MatchPredicate& matchClass) const;

        //- The names of objects with a matching class name
        //- that also have a matching object name.
        //  \note uses the class type() method
        template<class MatchPredicate1, class MatchPredicate2>
        wordList names
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The names of objects with a class satisfying \c isA\<Type\>.
        //
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type>
        wordList names() const;

        //- The names of objects with a class satisfying \c isA\<Type\>
        //- that also have a matching object name.
        //
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type, class MatchPredicate>
        wordList names(const MatchPredicate& matchName) const;


    // Summary of names (sorted)

        //- The sorted names of all objects
        wordList sortedNames() const;

        //- The sorted names of objects with the given class name.
        //  \note uses the class type() method
        wordList sortedNames(const char* clsName) const;

        //- The sorted names objects with a matching class name
        //  \note uses the class type() method
        template<class MatchPredicate>
        wordList sortedNames(const MatchPredicate& matchClass) const;

        //- The sorted names of objects with a matching class name
        //- that also have a matching object name.
        //  \note uses the class type() method
        template<class MatchPredicate1, class MatchPredicate2>
        wordList sortedNames
        (
            const MatchPredicate1& matchClass,
            const MatchPredicate2& matchName
        ) const;

        //- The sorted names of objects with a class satisfying \c isA\<Type\>
        //
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type>
        wordList sortedNames() const;

        //- The sorted names of objects with a class satisfying \c isA\<Type\>
        //- that also have a matching object name.
        //
        //  \note If \a Type is \c void, no isA check is used (always true).
        template<class Type, class MatchPredicate>
        wordList sortedNames(const MatchPredicate& matchName) const;


    // Lookup

        //- Lookup and return a const sub-objectRegistry.
        //
        //  \param forceCreate create it if it does not exist.
        //  \param recursive search parent registries.
        const objectRegistry& subRegistry
        (
            const word& name,
            const bool forceCreate = false,
            const bool recursive = false
        ) const;


        //- Return all objects with a class satisfying \c isA\<Type\>
        //
        //  \param strict use \c isType\<Type\> instead of \c isA\<Type\>
        template<class Type>
        HashTable<const Type*> lookupClass(const bool strict = false) const;

        //- Return all objects with a class satisfying \c isA\<Type\>
        //
        //  \param strict use \c isType\<Type\> instead of \c isA\<Type\>
        template<class Type>
        HashTable<Type*> lookupClass(const bool strict = false);

        //- Return const pointer to the regIOobject.
        //
        //  \param recursive search parent registries
        //
        //  \return nullptr if the object was not found.
        const regIOobject* cfindIOobject
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Can the regIOobject object be found (by name).
        //
        //  \param recursive search parent registries
        bool found(const word& name, const bool recursive = false) const;


        //- Is the named Type found?
        //
        //  \param recursive search parent registries
        template<class Type>
        bool foundObject
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Return const pointer to the object of the given Type.
        //
        //  \param recursive search parent registries
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class Type>
        const Type* cfindObject
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Return const pointer to the object of the given Type.
        //
        //  \param recursive search parent registries
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class Type>
        const Type* findObject
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Return non-const pointer to the object of the given Type.
        //
        //  \param recursive search parent registries
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class Type>
        Type* findObject
        (
            const word& name,
            const bool recursive = false
        );

        //- Return non-const pointer to the object of the given Type,
        //- using a const-cast to have it behave like a mutable.
        //  Exercise caution when using.
        //
        //  \param recursive search parent registries.
        //
        //  \return nullptr if the object was not found or had incorrect type.
        template<class Type>
        Type* getObjectPtr
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Lookup and return const reference to the object
        //- of the given Type. Fatal if not found or the wrong type.
        //
        //  \param recursive search parent registries.
        template<class Type>
        const Type& lookupObject
        (
            const word& name,
            const bool recursive = false
        ) const;

        //- Lookup and return non-const reference to the object
        //- of the given Type. Fatal if not found or the wrong type.
        //
        //  \param recursive search parent registries.
        template<class Type>
        Type& lookupObjectRef
        (
            const word& name,
            const bool recursive = false
        ) const;


    // Events

        //- Return new event number.
        label getEvent() const;


    // Edit

        //- Add a regIOobject to registry. A nullptr is ignored.
        bool checkIn(regIOobject* io) const;

        //- Add a regIOobject to registry
        bool checkIn(regIOobject& io) const;

        //- Remove a regIOobject from registry and free memory if the
        //- object is ownedByRegistry. A nullptr is ignored.
        bool checkOut(regIOobject* io) const;

        //- Remove a regIOobject from registry and free memory if the
        //- object is ownedByRegistry.
        bool checkOut(regIOobject& io) const;

        //- Remove a regIOobject by name from registry and free memory if the
        //- object is ownedByRegistry
        bool checkOut(const word& key) const;

        //- Clear all entries from the registry
        //  Performs a checkOut() for all objects that are ownedByRegistry
        void clear();

        //- Clear all entries from the registry and the table itself.
        void clearStorage();

        //- Erase an entry specified by the given iterator.
        //  Performs a checkOut() if the object was ownedByRegistry.
        //  \return True if the entry existed and was removed
        bool erase(const iterator& iter);

        //- Erase an entry specified by the given key
        //  Performs a checkOut() if the object was ownedByRegistry.
        //  \return True if the entry existed and was removed
        bool erase(const word& key);

        //- Remove entries given by the listed keys
        //  Performs a checkOut() for all objects that are ownedByRegistry.
        //  \return The number of items removed
        label erase(std::initializer_list<word> keys);

        //- Remove entries given by the listed keys
        //  Performs a checkOut() for all objects that are ownedByRegistry.
        //  \return The number of items removed
        label erase(const UList<word>& keys);

        //- Rename
        virtual void rename(const word& newName);


    // Reading

        //- Return true if any of the object's files have been modified
        virtual bool modified() const;

        //- Read the objects that have been modified
        void readModifiedObjects();

        //- Read object if modified
        virtual bool readIfModified();


    // Writing

        //- writeData function required by regIOobject but not used.
        //  For this class, write is used instead
        virtual bool writeData(Ostream&) const
        {
            NotImplemented;
            return false;
        }

        //- Write the objects using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool valid
        ) const;


    // Housekeeping

        //- Deprecated(2018-10) find object
        //  \deprecated(2018-10) - use findObject() method
        template<class Type>
        FOAM_DEPRECATED_FOR(2018-10, "findObject / cfindObject() methods")
        const Type* lookupObjectPtr
        (
            const word& name,
            bool recursive = false
        ) const
        {
            return this->cfindObject<Type>(name, recursive);
        }

        //- Deprecated(2018-10) get object pointer, ignoring constness
        //  \deprecated(2018-10) - use getObjectPtr() method
        template<class Type>
        FOAM_DEPRECATED_FOR(2018-10, "getObjectPtr() method")
        Type* lookupObjectRefPtr
        (
            const word& name,
            bool recursive = false
        ) const
        {
            return this->getObjectPtr<Type>(name, recursive);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "objectRegistryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
