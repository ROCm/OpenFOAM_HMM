/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
#include "dictionary.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Default is 2 : report source file name and line number if available
int Foam::messageStream::level(Foam::debug::infoSwitch("outputLevel", 2));

int Foam::messageStream::redirect(0);

// Default is 1 : report to Info
int Foam::infoDetailLevel(1);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::messageStream::messageStream
(
    const string& title,
    const errorSeverity severity,
    const int maxErrors
)
:
    title_(title),
    severity_(severity),
    maxErrors_(maxErrors),
    errorCount_(0)
{}


Foam::messageStream::messageStream(const dictionary& dict)
:
    title_(dict.get<string>("title")),
    severity_(FATAL),
    maxErrors_(0),
    errorCount_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::OSstream& Foam::messageStream::stream(OSstream* alternative)
{
    if (level)
    {
        // Serlal (master only) output?
        const bool serialOnly
        (
            (
                severity_ == INFO
             || severity_ == INFO_STDERR
             || severity_ == WARNING
            )
         || !UPstream::parRun()
        );

        if (serialOnly && (UPstream::parRun() && !UPstream::master()))
        {
            return Snull; // Non-serial, non-master: exit early
        }


        // Use stderr instead of stdout:
        // - requested via static <redirect> variable
        // - explicit:  INFO_STDERR
        // - inferred:  WARNING -> stderr when infoDetailLevel == 0
        const bool useStderr =
        (
            (redirect == 2)
         || (severity_ == INFO_STDERR)
         || (severity_ == WARNING && Foam::infoDetailLevel == 0)
        );

        OSstream* osptr;

        if (serialOnly)
        {
            // Use supplied alternative? Valid for serial only
            osptr = alternative;

            if (!osptr)
            {
                osptr = (useStderr ? &Serr : &Sout);
            }
        }
        else
        {
            // Non-serial
            osptr = (useStderr ? &Perr : &Pout);
        }

        if (!title_.empty())
        {
            (*osptr) << title_.c_str();
        }

        if (maxErrors_ && (++errorCount_ >= maxErrors_))
        {
            FatalErrorInFunction
                << "Too many errors..."
                << abort(FatalError);
        }

        return *osptr;
    }

    return Snull;
}


Foam::OSstream& Foam::messageStream::masterStream(const label communicator)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** messageStream with comm:" << communicator << endl;
        error::printStack(Pout);
    }

    if (communicator == UPstream::worldComm || UPstream::master(communicator))
    {
        return this->stream();
    }

    return Snull;
}


std::ostream& Foam::messageStream::stdStream()
{
    return this->stream().stdStream();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::OSstream& Foam::messageStream::operator()
(
    const string& functionName
)
{
    OSstream& os = this->stream();

    if (!functionName.empty())
    {
        os  << nl
            << "    From " << functionName.c_str() << nl;
    }

    return os;
}


Foam::OSstream& Foam::messageStream::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    OSstream& os = this->stream();

    os  << nl
        << "    From " << functionName << nl
        << "    in file " << sourceFileName
        << " at line " << sourceFileLineNumber << endl
        << "    ";

    return os;
}


Foam::OSstream& Foam::messageStream::operator()
(
    const string& functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    return operator()
    (
        functionName.c_str(),
        sourceFileName,
        sourceFileLineNumber
    );
}


Foam::OSstream& Foam::messageStream::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const string& ioFileName,
    const label ioStartLineNumber,
    const label ioEndLineNumber
)
{
    OSstream& os = this->stream();

    os  << nl
        << "    From " << functionName << nl
        << "    in file " << sourceFileName
        << " at line " << sourceFileLineNumber << nl
        << "    Reading " << ioFileName;

    if (ioStartLineNumber >= 0)
    {
        os  << " at line " << ioStartLineNumber;

        if (ioStartLineNumber < ioEndLineNumber)
        {
            os  << " to " << ioEndLineNumber;
        }
    }

    os << endl  << "    ";

    return os;
}


Foam::OSstream& Foam::messageStream::operator()
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
        -1  // No known endLineNumber
    );
}


Foam::OSstream& Foam::messageStream::operator()
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
        dict.relativeName(),
        dict.startLineNumber(),
        dict.endLineNumber()
    );
}


// * * * * * * * * * * * * * * * Global Variables  * * * * * * * * * * * * * //

Foam::messageStream Foam::Info("", Foam::messageStream::INFO);

Foam::messageStream Foam::InfoErr("", Foam::messageStream::INFO_STDERR);

Foam::messageStream Foam::Warning
(
    "--> FOAM Warning : ",
    Foam::messageStream::WARNING
);

Foam::messageStream Foam::SeriousError
(
    "--> FOAM Serious Error : ",
    Foam::messageStream::SERIOUS,
    100
);


// ************************************************************************* //
