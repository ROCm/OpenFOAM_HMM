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

\*---------------------------------------------------------------------------*/

// Foam header files.
#include "OSspecific.H"

// Project header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "LogManager.H"
#include "LogEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::LogManager* Foam::LogManager::pLogManager_ = NULL;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::LogManager::LogManager(const Foam::fileName& logFile)
:
    infoOpen_(false)
{
    static const char* functionName =
        "LogManager::LogManager(const Foam::fileName&)";

    try
    {
        // Check that there isn't already a global log manager object.
        if
        (
            FoamX::FoamXError::debug && Foam::LogManager::pLogManager_ == NULL
        )
        {
            Foam::LogManager::pLogManager_ = this;

            // Make sure the path exists.
            if (!dir(logFile.path()) && !mkDir(logFile.path()))
            {
                throw FoamX::FoamXError
                (
                    FoamXServer::E_FAIL,
                    "Log file directory '" +  logFile.path()
                  + "' could not be created.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Open the file stream.
            pFileStream_ = new Foam::OFstream(logFile);

            // Write XML header.
            (*pFileStream_) << "<?xml version=\"1.0\" ?>\n" << endl;
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::LogManager::~LogManager()
{
    if (Foam::LogManager::pLogManager_ == this)
    {
        Foam::LogManager::pLogManager_ = NULL;

        if (pFileStream_ != NULL)
        {
            delete pFileStream_;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::LogManager* Foam::LogManager::GlobalLogManager()
{
    return pLogManager_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::LogManager::startLogEntry(const LogEntry* pLogEntry)
{
    static const char* functionName =
        "LogManager::startLogEntry(const LogEntry*)";

    try
    {
        if (infoOpen_)
        {
            endInfoEntry();
        }

        // Push new entry onto stack.
        entryStack_.push(pLogEntry);

        pLogEntry->writeStart((*pFileStream_));
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::LogManager::endLogEntry()
{
    static const char* functionName = "LogManager::endLogEntry()";

    try
    {
        if (infoOpen_)
        {
            endInfoEntry();
        }

        // Pop the current entry off the stack.
        const LogEntry* pLogEntry = entryStack_.top();
        entryStack_.pop();

        pLogEntry->writeEnd((*pFileStream_));
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::LogManager::startInfoEntry()
{
    static const char* functionName =
        "LogManager::startInfoEntry()";

    try
    {
        if (infoOpen_)
        {
            endInfoEntry();
        }

        *pFileStream_ << "  <FoamXInfo>" << endl;
        infoOpen_ = true;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::LogManager::endInfoEntry()
{
    static const char* functionName =
        "LogManager::endInfoEntry()";

    try
    {
        *pFileStream_ << "  </FoamXInfo>" << endl;
        infoOpen_ = false;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::LogManager::operator Foam::OSstream&()
{
    static const char* functionName =
        "LogManager::operator Foam::OSstream&()";

    try
    {
        startInfoEntry();
        return *pFileStream_;
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
