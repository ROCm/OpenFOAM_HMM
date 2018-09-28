/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::TimePaths

Description
    Address the time paths without using the Time class.

SourceFiles
    TimePaths.C

\*---------------------------------------------------------------------------*/

#ifndef TimePaths_H
#define TimePaths_H

#include "fileName.H"
#include "instantList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class argList;

/*---------------------------------------------------------------------------*\
                           Class TimePaths Declaration
\*---------------------------------------------------------------------------*/


class TimePaths
{
    // Private data

        bool processorCase_;
        const fileName rootPath_;
        bool distributed_;
        fileName globalCaseName_;
        fileName case_;
        const word system_;
        const word constant_;


    // Private Member Functions

        //- Determine from case name whether it is a processor directory
        bool detectProcessorCase();


public:

    // Constructors

        //- Construct using characteristics given by the argList
        TimePaths
        (
            const argList& args,
            const word& systemName = "system",
            const word& constantName = "constant"
        );

        //- Construct given database name, rootPath and casePath
        TimePaths
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant"
        );


        //- Construct given database name, rootPath and casePath
        TimePaths
        (
            const bool processorCase,
            const fileName& rootPath,
            const bool distributed,
            const fileName& globalCaseName,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant"
        );


    // Member Functions

        //- Return true if this is a processor case
        bool processorCase() const
        {
            return processorCase_;
        }

        //- Return root path
        const fileName& rootPath() const
        {
            return rootPath_;
        }

        //- Return global case name
        const fileName& globalCaseName() const
        {
            return globalCaseName_;
        }

        //- Return case name
        const fileName& caseName() const
        {
            return case_;
        }

        //- Return case name
        fileName& caseName()
        {
            return case_;
        }

        //- Return system name
        const word& system() const
        {
            return system_;
        }

        //- Return the system name for the case, which is
        //- \c ../system() for parallel runs.
        fileName caseSystem() const;

        //- Return constant name
        const word& constant() const
        {
            return constant_;
        }

        //- Is case running with parallel distributed directories
        //  (i.e. not NFS mounted)
        bool distributed() const
        {
            return distributed_;
        }

        //- Return the constant name for the case, which is
        //- \c ../constant() for parallel runs.
        fileName caseConstant() const;

        //- Return path for the case
        fileName path() const
        {
            return rootPath()/caseName();
        }

        //- Return global path for the case
        fileName globalPath() const
        {
            return rootPath()/globalCaseName();
        }

        //- Return system path
        fileName systemPath() const
        {
            return path()/system();
        }

        //- Return constant path
        fileName constantPath() const
        {
            return path()/constant();
        }


    // Searching

        //- Search a given directory for valid time directories
        //  Forwards to the current fileHandler
        static instantList findTimes
        (
            const fileName& directory,
            const word& constantName = "constant"
        );

        //- Search instantList for the time index closest to the specified time
        static label findClosestTimeIndex
        (
            const instantList& timeDirs,
            const scalar t,
            const word& constantName = "constant"
        );

        //- Search the case for valid time directories
        instantList times() const;

        //- Search the case for the time closest to the given time
        instant findClosestTime(const scalar t) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
