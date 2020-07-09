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
    Foam::profilingInformation

Description
    Code profiling information in terms of time spent, number of calls etc.

SourceFiles
    profilingInformation.C

\*---------------------------------------------------------------------------*/

#ifndef profilingInformation_H
#define profilingInformation_H

#include "label.H"
#include "scalar.H"
#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class profiling;

/*---------------------------------------------------------------------------*\
                    Class profilingInformation Declaration
\*---------------------------------------------------------------------------*/

class profilingInformation
{
    // Private Data

        //- Unique id to identify it
        const label id_;

        //- What this timer does
        const string description_;

        //- Pointer to the parent object (or self for top-level)
        profilingInformation* parent_;

        //- Number of times this was called
        long calls_;

        //- Total time spent
        scalar totalTime_;

        //- Time spent in children
        scalar childTime_;

        //- Max memory usage on call.
        //  Only valid when the calling profiling has memInfo active.
        mutable int maxMem_;

        //- Is this information active or passive (ie, on the stack)?
        mutable bool active_;


protected:

    // Friendship

        //- Allow creation of master-element, setActive etc.
        friend class profiling;


    // Constructors

        //- Default construct - only the master-element
        profilingInformation();


    // Protected Member Functions

        //- Mark as being active or passive)
        void setActive(bool state) const;

        //- No copy construct
        profilingInformation(const profilingInformation&) = delete;

        //- No copy assignment
        void operator=(const profilingInformation&) = delete;

public:

    // Constructors

        //- Construct from components
        profilingInformation
        (
            profilingInformation* parent,
            const string& descr,
            const label id
        );


    //- Destructor
    ~profilingInformation() = default;


    // Member Functions

    // Access

        label id() const
        {
            return id_;
        }

        const string& description() const
        {
            return description_;
        }

        profilingInformation& parent() const
        {
            return *parent_;
        }

        label calls() const
        {
            return calls_;
        }

        scalar totalTime() const
        {
            return totalTime_;
        }

        scalar childTime() const
        {
            return childTime_;
        }

        int maxMem() const
        {
            return maxMem_;
        }

        bool active() const
        {
            return active_;
        }


    // Edit

        //- Update it with a new timing information
        void update(const scalar elapsedTime);


    // Write

        //- Write the profiling times, optionally with additional values
        //  Uses dictionary format.
        Ostream& write
        (
            Ostream& os,
            const bool offset = false,
            const scalar elapsedTime = 0,
            const scalar childTime = 0
        ) const;
};


// Global Operators

Ostream& operator<<(Ostream& os, const profilingInformation& info);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
