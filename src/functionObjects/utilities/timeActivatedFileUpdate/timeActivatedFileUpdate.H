/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::functionObjects::timeActivatedFileUpdate

Group
    grpUtilitiesFunctionObjects

Description
    Performs a file copy/replacement once a specified time has been reached.

Usage
    Example usage to update the fvSolution dictionary at various times
    throughout the calculation:
    \verbatim
    fileUpdate1
    {
        type              timeActivatedFileUpdate;
        libs              ("libutilityFunctionObjects.so");
        writeControl      timeStep;
        writeInterval     1;
        fileToUpdate      "$FOAM_CASE/system/fvSolution";
        timeVsFile
        (
            (-1 "$FOAM_CASE/system/fvSolution.0")
            (0.10 "$FOAM_CASE/system/fvSolution.10")
            (0.20 "$FOAM_CASE/system/fvSolution.20")
            (0.35 "$FOAM_CASE/system/fvSolution.35")
        );
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: timeActivatedFileUpdate | yes |
        fileToUpdate | Name of file to update  | yes         |
        timeVsFile   | List of time vs file    | yes         |
    \endtable


SourceFiles
    timeActivatedFileUpdate.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_timeActivatedFileUpdate_H
#define functionObjects_timeActivatedFileUpdate_H

#include "functionObject.H"
#include "Tuple2.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class timeActivatedFileUpdate Declaration
\*---------------------------------------------------------------------------*/

class timeActivatedFileUpdate
:
    public functionObject
{
    // Private data

        //- Reference to Time
        const Time& time_;

        //- Name of file to update
        fileName fileToUpdate_;

        //- List of times vs filenames
        List<Tuple2<scalar, fileName>> timeVsFile_;

        //- Index of last file copied
        label lastIndex_;

        //- Has anything been copied?
        bool modified_;


    // Private Member Functions

        //- Update file
        void updateFile();

        //- Disallow default bitwise copy construct
        timeActivatedFileUpdate(const timeActivatedFileUpdate&);

        //- Disallow default bitwise assignment
        void operator=(const timeActivatedFileUpdate&);


public:

    //- Runtime type information
    TypeName("timeActivatedFileUpdate");


    // Constructors

        //- Construct from Time and dictionary
        timeActivatedFileUpdate
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~timeActivatedFileUpdate();


    // Member Functions

        //- Read the timeActivatedFileUpdate data
        virtual bool read(const dictionary&);

        //- Execute file updates
        virtual bool execute();

        //- Do nothing
        virtual bool write();

        //- Did any file get changed during execution?
        virtual bool filesModified() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
