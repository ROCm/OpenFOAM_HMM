/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

// Forward Declarations
class argList;

/*---------------------------------------------------------------------------*\
                           Class TimePaths Declaration
\*---------------------------------------------------------------------------*/


class TimePaths
{
    // Private Data

        bool processorCase_;
        bool distributed_;

        const fileName rootPath_;
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
        explicit TimePaths
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

        //- True if case running with parallel distributed directories
        //- (ie. not NFS mounted)
        inline bool distributed() const noexcept;

        //- Return true if this is a processor case
        inline bool processorCase() const noexcept;

        //- Set processor case status. Return old one
        inline bool processorCase(const bool isProcessorCase) noexcept;

        //- Return root path
        inline const fileName& rootPath() const;

        //- Return global case name
        inline const fileName& globalCaseName() const;

        //- Return case name
        inline const fileName& caseName() const;

        //- The case name for modification (use with caution)
        inline fileName& caseName();

        //- Return path for the case
        inline fileName path() const;

        //- Return global path for the case
        inline fileName globalPath() const;

        //- Return the input relative to the globalPath by stripping off
        //- a leading value of the globalPath
        //
        //  \param input the directory or filename to make case-relative
        //  \param caseTag replace globalPath with \<case\> for later
        //      use with expand(), or prefix \<case\> if the file name was
        //      not an absolute location
        inline fileName relativePath
        (
            const fileName& input,
            const bool caseTag = false
        ) const;


        //- Return constant name
        inline const word& constant() const;

        //- Return system name
        inline const word& system() const;

        //- Return the constant name for the case, which is
        //- \c ../constant() for parallel runs.
        inline fileName caseConstant() const;

        //- Return the system name for the case, which is
        //- \c ../system() for parallel runs.
        inline fileName caseSystem() const;

        //- Return constant path
        inline fileName constantPath() const;

        //- Return system path
        inline fileName systemPath() const;


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

#include "TimePathsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
