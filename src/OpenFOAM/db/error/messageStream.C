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

Foam::OSstream& Foam::messageStream::masterStream(const label communicator)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** messageStream with comm:" << communicator << endl;
        error::printStack(Pout);
    }

    if (communicator == UPstream::worldComm || UPstream::master(communicator))
    {
        return operator()();
    }

    return Snull;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::OSstream& Foam::messageStream::operator()
(
    const string& functionName
)
{
    OSstream& os = operator OSstream&();

    os  << nl
        << "    From " << functionName.c_str() << nl;

    return os;
}


Foam::OSstream& Foam::messageStream::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    OSstream& os = operator OSstream&();

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
    OSstream& os = operator OSstream&();

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
        -1
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
        dict.name(),
        dict.startLineNumber(),
        dict.endLineNumber()
    );
}


Foam::messageStream::operator Foam::OSstream&()
{
    if (level)
    {
        const bool collect =
        (
            severity_ == INFO
         || severity_ == WARNING
         || severity_ == INFO_STDERR
        );

        // Could add guard with parRun
        if (collect && !Pstream::master())
        {
            return Snull;
        }

        // Use stderr instead of stdout
        // - INFO_STDERR
        // - WARNING when infoDetailLevel == 0
        const bool useStderr =
        (
            (severity_ == INFO_STDERR)
         || (severity_ == WARNING && Foam::infoDetailLevel == 0)
        );

        OSstream& os =
        (
            (collect || !Pstream::parRun())
          ? (useStderr ? Serr : Sout)
          : (useStderr ? Perr : Pout)
        );


        if (!title().empty())
        {
            os << title().c_str();
        }

        if (maxErrors_ && (++errorCount_ >= maxErrors_))
        {
            FatalErrorInFunction
                << "Too many errors"
                << abort(FatalError);
        }

        return os;
    }

    return Snull;
}


// * * * * * * * * * * * * * * * Global Variables  * * * * * * * * * * * * * //

Foam::messageStream Foam::Info("", messageStream::INFO);

Foam::messageStream Foam::InfoErr("", messageStream::INFO_STDERR);

Foam::messageStream Foam::Warning
(
    "--> FOAM Warning : ",
    messageStream::WARNING
);

Foam::messageStream Foam::SeriousError
(
    "--> FOAM Serious Error : ",
    messageStream::SERIOUS,
    100
);


// ************************************************************************* //
