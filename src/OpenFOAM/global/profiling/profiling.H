/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009-2016 Bernhard Gschaider
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

    This is typically activated from within system/controlDict as follows
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
#include "DynamicList.H"
#include "PtrDynList.H"
#include "Time.H"
#include "clockTime.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Ostream;
class cpuInfo;
class memInfo;
class profilingInformation;
class profilingSysInfo;
class profilingTrigger;

/*---------------------------------------------------------------------------*\
                          Class profiling Declaration
\*---------------------------------------------------------------------------*/

class profiling
:
    public IOdictionary
{
public:

    // Public Typedefs

        typedef profilingInformation Information;
        typedef profilingTrigger Trigger;


    // Static Data Members

        //- Flag if profiling is allowed
        static int allowed;

private:

    // Private Typedefs

        typedef profilingSysInfo sysInfo;


    // Private Static Data Members

        //- Only one global object is possible
        static std::unique_ptr<profiling> singleton_;


    // Private Data Members

        //- The owner of the profiling
        const Time& owner_;

        //- Storage of profiling information (memory management)
        PtrDynList<Information> pool_;

        //- Parent/child relationships for lookup purposes
        DynamicList<DynamicList<Information*,16>> children_;

        //- LIFO stack of profiling information
        DynamicList<Information*> stack_;

        //- LIFO stack of clock values
        DynamicList<clockValue> times_;

        //- General system information (optional)
        std::unique_ptr<sysInfo> sysInfo_;

        //- CPU-Information (optional)
        std::unique_ptr<cpuInfo> cpuInfo_;

        //- MEM-Information (optional)
        std::unique_ptr<memInfo> memInfo_;


    // Private Member Functions

        //- No copy construct
        profiling(const profiling&) = delete;

        //- No copy assignment
        void operator=(const profiling&) = delete;


protected:

    // Friendship

        friend class profilingTrigger;
        friend class Time;


    // Constructors

        //- Construct IO object, everything enabled
        profiling
        (
            const IOobject& io,
            const Time& owner,
            const bool allEnabled = true
        );

        //- Construct IO object with finer control over behaviour
        profiling
        (
            const dictionary& dict,
            const IOobject& io,
            const Time& owner
        );


    // Protected Member Functions

        //- Clear all profiling and restart with new profiling
        //  \return pointer to stored information element
        Information* create();

        //- Get or create named profiling information element with the
        //- specified parent.
        //  \return pointer to stored information element
        Information* create
        (
            Information* parent,
            const string& descr
        );

        //- Add to stack of active information and begin timer datum
        void beginTimer(Information* info);

        //- Remove from stack of active information and update elapsed time
        //  \return pointer to profiling information element (for reference)
        Information* endTimer();


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
        //  Returns nullptr if profiling has not been initialized
        static profilingInformation* New(const string& descr);

        //- Remove the information from the top of the stack
        static void unstack(const profilingInformation* info);


public:

    //- Destructor. Top-level clears the singleton.
    ~profiling();


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
        label size() const noexcept;

        //- writeData member function required by regIOobject
        virtual bool writeData(Ostream& os) const;

        //- Write as uncompressed ASCII
        virtual bool writeObject
        (
            IOstreamOption /*ignore*/,
            const bool valid
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
