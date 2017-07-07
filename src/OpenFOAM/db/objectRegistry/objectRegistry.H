/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
        bool parentNotTime() const;

        //- Templated implementation for classes()
        template<class UnaryMatchPredicate>
        static HashTable<wordHashSet> classesImpl
        (
            const objectRegistry& list,
            const UnaryMatchPredicate& matcher
        );

        //- Templated implementation for names()
        template<class Type, class UnaryMatchPredicate>
        static wordList namesImpl
        (
            const objectRegistry& list,
            const UnaryMatchPredicate& matcher,
            const bool doSort
        );


        //- Disallow Copy constructor
        objectRegistry(const objectRegistry&) = delete;

        //- Disallow default bitwise copy construct and assignment
        void operator=(const objectRegistry&) = delete;


public:

    //- Declare type name for this IOobject
    TypeName("objectRegistry");


    // Constructors

        //- Construct the time objectRegistry given an initial estimate
        //  for the number of entries
        explicit objectRegistry
        (
            const Time& db,
            const label nIoObjects = 128
        );

        //- Construct a sub-registry given an IObject to describe the registry
        //  and an initial estimate for the number of entries
        explicit objectRegistry
        (
            const IOobject& io,
            const label nIoObjects = 128
        );


    //- Destructor
    virtual ~objectRegistry();


    // Member functions

        // Access

            //- Return time
            const Time& time() const
            {
                return time_;
            }

            //- Return the parent objectRegistry
            const objectRegistry& parent() const
            {
                return parent_;
            }

            //- Local directory path of this objectRegistry relative to the time
            virtual const fileName& dbDir() const
            {
                return dbDir_;
            }


        // Summary of classes

            //- A summary hash of classes used and their associated object names.
            //  Behaviour and usage as per IOobjectList::classes
            HashTable<wordHashSet> classes() const;

            //- A summary hash of classes used and their associated object names
            //  restricted to objects with names that satisfy the input matcher
            HashTable<wordHashSet> classes(const wordRe& matcher) const;

            //- A summary hash of classes used and their associated object names
            //  restricted to objects with names that satisfy the input matcher
            HashTable<wordHashSet> classes(const wordRes& matcher) const;


        // Summary of names

            //- A list of names of the objects
            wordList names() const;

            //- The names of objects with the given class name
            wordList names(const word& clsName) const;

            //- The names of objects with the given type
            template<class Type>
            wordList names() const;

            //- The names of objects with the given type that also
            //  have a name satisfying the input matcher
            template<class Type>
            wordList names(const wordRe& matcher) const;

            //- The names of objects with the given type that also
            //  have a name satisfying the input matcher
            template<class Type>
            wordList names(const wordRes& matcher) const;


        // Summary of names (sorted)

            //- A sorted list of names of the objects
            wordList sortedNames() const;

            //- The sorted names of objects with the given class name
            wordList sortedNames(const word& clsName) const;

            //- The sorted names of objects with the given type
            template<class Type>
            wordList sortedNames() const;

            //- The sorted names of objects with the given type that also
            //  have a name satisfying the input matcher
            template<class Type>
            wordList sortedNames(const wordRe& matcher) const;

            //- The sorted names of objects with the given type that also
            //  have a name satisfying the input matcher
            template<class Type>
            wordList sortedNames(const wordRes& matcher) const;


        // Lookup

            //- Lookup and return a const sub-objectRegistry.
            //  Optionally create it if it does not exist.
            //  If recursive, search parent registries.
            const objectRegistry& subRegistry
            (
                const word& name,
                const bool forceCreate = false,
                const bool recursive = false
            ) const;


            //- Lookup and return all objects of the given Type
            template<class Type>
            HashTable<const Type*> lookupClass(const bool strict = false) const;

            //- Lookup and return all objects of the given Type
            template<class Type>
            HashTable<Type*> lookupClass(const bool strict = false);

            //- Is the named Type found?
            //  If recursive, search parent registries.
            template<class Type>
            bool foundObject
            (
                const word& name,
                const bool recursive = false
            ) const;

            //- Lookup and return the object of the given Type.
            //  If recursive, search parent registries.
            template<class Type>
            const Type& lookupObject
            (
                const word& name,
                const bool recursive = false
            ) const;

            //- Lookup and return the object of the given Type.
            //  If recursive, search parent registries.
            template<class Type>
            Type& lookupObjectRef
            (
                const word& name,
                const bool recursive = false
            ) const;

            //- Lookup and return pointer to the object of the given Type,
            //  otherwise nullptr if the object was not found,
            //  or had the incorrect type.
            //  If recursive, search parent registries.
            template<class Type>
            const Type* lookupObjectPtr
            (
                const word& name,
                const bool recursive = false
            ) const;


            //- Lookup and return non-const pointer to the object
            //  of the given Type,
            //  otherwise nullptr if the object was not found,
            //  or had the incorrect type.
            //  If recursive, search parent registries.
            template<class Type>
            Type* lookupObjectRefPtr
            (
                const word& name,
                const bool recursive = false
            ) const;


        // Events

            //- Return new event number.
            label getEvent() const;


        // Edit

            //- Rename
            virtual void rename(const word& newName);

            //- Add an regIOobject to registry
            bool checkIn(regIOobject& io) const;

            //- Remove an regIOobject from registry
            bool checkOut(regIOobject& io) const;


        // Reading

            //- Return true if any of the object's files have been modified
            virtual bool modified() const;

            //- Read the objects that have been modified
            void readModifiedObjects();

            //- Read object if modified
            virtual bool readIfModified();


        // Writing

            //- writeData function required by regIOobject but not used
            //  for this class, write is used instead
            virtual bool writeData(Ostream&) const
            {
                NotImplemented;

                return false;
            }

            //- Write the objects
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool valid
            ) const;
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
