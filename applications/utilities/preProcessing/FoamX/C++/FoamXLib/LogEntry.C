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

// Project header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "LogManager.H"
#include "LogEntry.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::LogEntry::LogEntry
(
    const char* methodName,
    const char* fileName,
    const int lineNo
)
:
    logManagerPtr_(NULL),
    lineNo_(-1)
{
    if (!Pstream::parRun() || Pstream::master())
    {
        // Get global log manager, if available.
        logManagerPtr_ = LogManager::GlobalLogManager();

        if (logManagerPtr_)
        {
            methodName_ = methodName;
            //methodName_.replaceAll("&", "&amp");
            methodName_.replaceAll("&", "");

            fileName_   = fileName;
            lineNo_     = lineNo;
            logManagerPtr_->startLogEntry(this);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::LogEntry::~LogEntry()
{
    if (logManagerPtr_)
    {
        logManagerPtr_->endLogEntry();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::LogEntry::writeStart(Ostream& os) const
{
    static const char* functionName =
        "LogEntry::writeStart(Ostream&) const";

    try
    {
        // Write opening tag.
        os << "<FunctionCall>" << endl;
        os << "  <MethodName>" << methodName_ << "</MethodName>" << endl;

        if (fileName_.size() > 0)
        {
            os << "  <FileName>" << fileName_ << "</FileName>" << endl;
        }

        if (lineNo_ > -1)
        {
            os << "  <LineNo>" << lineNo_ << "</LineNo>" << endl;
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::LogEntry::writeEnd(Ostream& os) const
{
    static const char* functionName =
        "LogEntry::writeEnd(Ostream&) const";

    try
    {
        // Write closing tag.
        os << "</FunctionCall>" << endl;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::LogEntry::operator Foam::OSstream&()
{
    if (logManagerPtr_)
    {
        return *logManagerPtr_;
    }
    else
    {
        return Snull;
    }
}


// ************************************************************************* //
