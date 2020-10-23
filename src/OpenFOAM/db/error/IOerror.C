/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "StringStream.H"
#include "fileName.H"
#include "dictionary.H"
#include "JobInfo.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOerror::IOerror(const string& title)
:
    error(title),
    ioFileName_("unknown"),
    ioStartLineNumber_(-1),
    ioEndLineNumber_(-1)
{}


Foam::IOerror::IOerror(const dictionary& errDict)
:
    error(errDict),
    ioFileName_(errDict.get<string>("ioFileName")),
    ioStartLineNumber_(errDict.get<label>("ioStartLineNumber")),
    ioEndLineNumber_(errDict.get<label>("ioEndLineNumber"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IOerror::~IOerror() noexcept
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::OSstream& Foam::IOerror::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const string& ioFileName,
    const label ioStartLineNumber,
    const label ioEndLineNumber
)
{
    error::operator()(functionName, sourceFileName, sourceFileLineNumber);
    ioFileName_ = ioFileName;
    ioStartLineNumber_ = ioStartLineNumber;
    ioEndLineNumber_ = ioEndLineNumber;

    return operator OSstream&();
}


Foam::OSstream& Foam::IOerror::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const IOstream& ioStream
)
{
    return operator()
    (
        functionName,
        sourceFileName,
        sourceFileLineNumber,
        ioStream.name(),
        ioStream.lineNumber(),
        -1
    );
}


Foam::OSstream& Foam::IOerror::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const dictionary& dict
)
{
    return operator()
    (
        functionName,
        sourceFileName,
        sourceFileLineNumber,
        dict.name(),
        dict.startLineNumber(),
        dict.endLineNumber()
    );
}


void Foam::IOerror::SafeFatalIOError
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const IOstream& ioStream,
    const string& msg
)
{
    if (JobInfo::constructed)
    {
        FatalIOError
        (
            functionName,
            sourceFileName,
            sourceFileLineNumber,
            ioStream
        )   << msg << Foam::exit(FatalIOError);
    }
    else
    {
        std::cerr
            << nl
            << "--> FOAM FATAL IO ERROR:" << nl
            << msg << nl
            << "file: " << ioStream.name()
            << " at line " << ioStream.lineNumber() << '.' << nl << nl
            << "    From " << functionName << nl
            << "    in file " << sourceFileName
            << " at line " << sourceFileLineNumber << '.' << std::endl;
        std::exit(1);
    }
}


Foam::IOerror::operator Foam::dictionary() const
{
    dictionary errDict(error::operator dictionary());

    errDict.remove("type");
    errDict.add("type", word("Foam::IOerror"));

    errDict.add("ioFileName", ioFileName());
    errDict.add("ioStartLineNumber", ioStartLineNumber());
    errDict.add("ioEndLineNumber", ioEndLineNumber());

    return errDict;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::IOerror::exitOrAbort(const int errNo, const bool isAbort)
{
    if (!throwing_ && JobInfo::constructed)
    {
        jobInfo.add("FatalIOError", operator dictionary());
        if (isAbort || error::useAbort())
        {
            jobInfo.abort();
        }
        else
        {
            jobInfo.exit();
        }
    }

    if (throwing_ && !isAbort)
    {
        // Make a copy of the error to throw
        IOerror errorException(*this);

        // Reset the message buffer for the next error message
        messageStreamPtr_->reset();

        throw errorException;
    }
    else if (error::useAbort())
    {
        Perr<< nl << *this << nl
            << "\nFOAM aborting (FOAM_ABORT set)\n" << endl;
        error::printStack(Perr);
        std::abort();
    }
    else if (Pstream::parRun())
    {
        if (isAbort)
        {
            Perr<< nl << *this << nl
                << "\nFOAM parallel run aborting\n" << endl;
            error::printStack(Perr);
            Pstream::abort();
        }
        else
        {
            Perr<< nl << *this << nl
                << "\nFOAM parallel run exiting\n" << endl;
            Pstream::exit(errNo);
        }
    }
    else
    {
        if (isAbort)
        {
            Perr<< nl << *this << nl
                << "\nFOAM aborting\n" << endl;
            error::printStack(Perr);

            #ifdef _WIN32
            std::exit(1);  // Prefer exit() to avoid unnecessary warnings
            #else
            std::abort();
            #endif
        }
        else
        {
            Perr<< nl << *this << nl
                << "\nFOAM exiting\n" << endl;
            std::exit(errNo);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IOerror::exit(const int)
{
    exitOrAbort(1, false);
}


void Foam::IOerror::abort()
{
    exitOrAbort(1, true);
}


void Foam::IOerror::write(Ostream& os, const bool includeTitle) const
{
    if (os.bad())
    {
        return;
    }

    os  << nl;
    if (includeTitle && !title().empty())
    {
        os  << title().c_str()
            << "(openfoam-" << foamVersion::api;

        if (foamVersion::patched())
        {
            // Patch-level, when defined
            os  << " patch=" << foamVersion::patch.c_str();
        }
        os  << ')' << nl;
    }
    os  << message().c_str();


    if (!ioFileName().empty())
    {
        os  << nl << nl
            << "file: " << ioFileName().c_str();

        if (ioStartLineNumber() >= 0)
        {
            os  << " at line " << ioStartLineNumber();
            if (ioStartLineNumber() < ioEndLineNumber())
            {
                os  << " to " << ioEndLineNumber();
            }
            os  << '.';
        }
    }


    const label lineNo = sourceFileLineNumber();

    if (IOerror::level >= 2 && lineNo && !functionName().empty())
    {
        os  << nl << nl
            << "    From " << functionName().c_str() << nl;

        if (!sourceFileName().empty())
        {
            os << "    in file " << sourceFileName().c_str();

            if (lineNo > 0)
            {
                os  << " at line " << lineNo << '.';
            }
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const IOerror& err)
{
    err.write(os);

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error definitions

Foam::IOerror Foam::FatalIOError("--> FOAM FATAL IO ERROR: ");


// ************************************************************************* //
