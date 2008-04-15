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
    LogManager.C

\*---------------------------------------------------------------------------*/

#ifndef LogManager_H
#define LogManager_H

// Foam header files.
#include "LIFOStack.H"
#include "fileName.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes.
class LogEntry;

/*---------------------------------------------------------------------------*\
                           Class LogManager Declaration
\*---------------------------------------------------------------------------*/

class LogManager
{

public:

    enum LogFilter
    {
        LOG_DEBUG,
        LOG_DIAGNOSTIC,
        LOG_WARNINGS,
        LOG_ERRORS 
    };


private:

    // Private data

        static LogManager* pLogManager_;
        LogFilter filter_;
        OFstream* pFileStream_;
        LIFOStack<const LogEntry*> entryStack_;
        bool infoOpen_;


public:

    // Constructors

        LogManager(const fileName& logFile);


    // Destructor

        ~LogManager();


    // Member Functions

        //- Access to the global log manager.
        static LogManager* GlobalLogManager();


        // Log entry methods.

            void startLogEntry(const LogEntry* pLogEntry);
            void endLogEntry();

            void startInfoEntry();
            void endInfoEntry();


        //- Convert to Ostream
        //  Prints basic message and then returns Ostream for further info.
        operator OSstream&();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
