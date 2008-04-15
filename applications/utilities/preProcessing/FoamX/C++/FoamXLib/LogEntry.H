/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::LogManager

Description

SourceFiles
    LogEntry.C

\*---------------------------------------------------------------------------*/

#ifndef LogEntry_H
#define LogEntry_H

// Foam header files.
#include "label.H"
#include "word.H"
#include "fileName.H"
#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class LogManager;

/*---------------------------------------------------------------------------*\
                           Class LogEntry Declaration
\*---------------------------------------------------------------------------*/

class LogEntry
{    
    // Private data

        LogManager* logManagerPtr_;
        string     methodName_;
        fileName   fileName_;
        label      lineNo_;

public:

    // Constructors

        LogEntry
        (
            const char* methodName,
            const char* fileName,
            const int lineNo
        );


    // Destructor

        ~LogEntry();


    // Member Functions

        void writeStart(Ostream&) const;

        void writeEnd(Ostream&) const;

        //- Convert to Ostream
        //  Prints basic message and then returns Ostream for further info.
        operator OSstream&();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
