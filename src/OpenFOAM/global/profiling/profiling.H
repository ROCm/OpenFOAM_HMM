/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2016 Bernhard Gschaider
     \\/     M anipulation  | Copyright (C) 2016-2107 OpenCFD Ltd.
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
    Foam::profiling

Description
    Code profiling.

    This is typically activated from within the system/controlDict as follows
    (defaults shown):
    \code
        profiling
        {
            active      true;
            cpuInfo     false;
            memInfo     false;
            sysInfo     false;
        }
    \endcode
    or simply using all defaults:
    \code
        profiling
        {}
    \endcode

SourceFiles
    profiling.C

\*---------------------------------------------------------------------------*/

#ifndef profiling_H
#define profiling_H

#include "profilingTrigger.H"
#include "IOdictionary.H"
#include "HashPtrTable.H"
#include "Tuple2.H"
#include "LIFOStack.H"
#include "Map.H"
#include "Time.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;
class cpuInfo;
class memInfo;
class profilingSysInfo;

/*---------------------------------------------------------------------------*\
                          Class profiling Declaration
\*---------------------------------------------------------------------------*/

class profiling
:
    public IOdictionary
{
public:

    // Public typedefs

        typedef profilingInformation Information;
        typedef profilingTrigger Trigger;

    // Static data members

        //- Flag if profiling is allowed
        static int allowed;

private:

    // Private classes, typedefs

        typedef profilingSysInfo sysInfo;

        //- Profile information lookup is qualified by parent id
        typedef Tuple2<string, label> Key;

        //- Hashing for information lookup
        class HashKey
        :
            public Hash<Key>
        {
        public:

            HashKey()
            {}

            //- Hash qualified by the parent id to avoid collisions
            unsigned operator()(const Key& key) const
            {
                return
                (
                    Hash<string>()(key.first())
                  + Hash<label>()(key.second())
                );
            }
        };

        typedef HashPtrTable<Information, Key, HashKey> StorageContainer;
        typedef LIFOStack<Information*> StackContainer;


    // Private Static Data Members

        //- Only one global pool object is possible
        static profiling* pool_;


    // Private Data Members

        //- The owner of the profiling
        const Time& owner_;

        //- A global timer for the profiling
        clockTime clockTime_;

        //- Storage of profiling information
        StorageContainer hash_;

        //- Local stack of profiling information
        StackContainer stack_;

        //- Note the timers (by Id) for the correct stack-output
        Map<clockTime*> timers_;

        //- General system information (optional)
        sysInfo* sysInfo_;

        //- CPU-Information (optional)
        cpuInfo* cpuInfo_;

        //- MEM-Information (optional)
        memInfo* memInfo_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        profiling(const profiling&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const profiling&) = delete;


protected:

    // Friendship

        friend class profilingTrigger;
        friend class Time;


    // Constructors

        //- Construct IO object, everything enabled
        profiling(const IOobject& io, const Time& owner);

        //- Construct IO object with finer control over behaviour
        profiling
        (
            const dictionary& dict,
            const IOobject& io,
            const Time& owner
        );


    //- Destructor
    ~profiling();


    // Protected Member Functions

        //- Find named profiling information element with specified parent.
        //  Return nullptr on failure.
        profilingInformation* find(const string& descr, const label parentId);

        //- Add to hashed storage,
        //  \return pointer to newly stored element for chaining
        profilingInformation* store(profilingInformation* info);

        //- Add to stack and set timer lookup (based on Id)
        void push(profilingInformation* info, clockTime& timer);

        //- Remove from stack and remove timer lookup (based on Id).
        //  \return pointer to profiling information element
        profilingInformation* pop();


    // Static control elements

        //- Singleton to initialize profiling pool, everything enabled
        static void initialize
        (
            const IOobject& ioObj,
            const Time& owner
        );

        //- Singleton to initialize profiling pool with finer control
        static void initialize
        (
            const dictionary& dict,
            const IOobject& ioObj,
            const Time& owner
        );

        //- Stop profiling, cleanup pool if possible
        static void stop(const Time& owner);

        //- Existing or new element on pool, add to stack.
        //  Returns null if profiling has not been initialized
        static profilingInformation* New
        (
            const string& descr,
            clockTime& timer
        );

        //- Remove the information from the top of the stack
        static void unstack(const profilingInformation* info);

public:

    // Static Member Functions

        //- True if profiling is allowed and is active
        static bool active();

        //- Disallow profiling by forcing the InfoSwitch off.
        static void disable();

        //- Print profiling information to specified output
        //  Forwards to writeData member of top-level object
        static bool print(Ostream& os);

        //- Write profiling information now
        static bool writeNow();


    // Member Functions

        //- The owner of the profiling
        const Time& owner() const;

        //- The size of the current stack
        label size() const;

        //- writeData member function required by regIOobject
        virtual bool writeData(Ostream& os) const;

        //- Write as uncompressed ASCII, using given format
        virtual bool writeObject
        (
            IOstream::streamFormat ignoreAlwaysASCII,
            IOstream::versionNumber ver,
            IOstream::compressionType ignoreAlwaysUncompressed,
            const bool valid
        ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
