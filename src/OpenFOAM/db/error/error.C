/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
#include "foamVersion.H"
#include "OSspecific.H"
#include "Switch.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::error::warnAboutAge(const char* what, const int version)
{
    if (version <= 0)
    {
        // No warning for 0 (unversioned) or -ve values (silent versioning)
    }
    else if (version < 1000)
    {
        // Warning for things that predate the YYMM versioning
        // (eg, 240 for version 2.4)
        std::cerr
            << "    This " << what << " is very old.\n"
            << std::endl;
    }
    else if (version < foamVersion::api)
    {
        const int months =
        (
            // YYMM -> months
            (12 * (foamVersion::api/100) + (foamVersion::api % 100))
          - (12 * (version/100)  + (version % 100))
        );

        std::cerr
            << "    This " << what << " is " << months << " months old.\n"
            << std::endl;
    }
    // No warning for (foamVersion::api < version).
    // We use this to denote future expiry dates of transition features.
}


bool Foam::error::useAbort()
{
    // FOAM_ABORT env set and contains bool-type value
    return static_cast<bool>(Switch::find(Foam::getEnv("FOAM_ABORT")));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::error::error(const string& title)
:
    std::exception(),
    messageStream(title, messageStream::FATAL),
    functionName_("unknown"),
    sourceFileName_("unknown"),
    sourceFileLineNumber_(0),
    throwing_(false),
    messageStreamPtr_(new OStringStream())
{}


Foam::error::error(const dictionary& errDict)
:
    std::exception(),
    messageStream(errDict),
    functionName_(errDict.get<string>("functionName")),
    sourceFileName_(errDict.get<string>("sourceFileName")),
    sourceFileLineNumber_(errDict.get<label>("sourceFileLineNumber")),
    throwing_(false),
    messageStreamPtr_(new OStringStream())
{}


Foam::error::error(const error& err)
:
    std::exception(),
    messageStream(err),
    functionName_(err.functionName_),
    sourceFileName_(err.sourceFileName_),
    sourceFileLineNumber_(err.sourceFileLineNumber_),
    throwing_(err.throwing_),
    messageStreamPtr_(new OStringStream(*err.messageStreamPtr_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::error::~error() noexcept
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::OSstream& Foam::error::operator()
(
    const string& functionName
)
{
    functionName_ = functionName;
    sourceFileName_.clear();
    sourceFileLineNumber_ = -1;

    return operator OSstream&();
}


Foam::OSstream& Foam::error::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    functionName_ = functionName;
    sourceFileName_ = sourceFileName;
    sourceFileLineNumber_ = sourceFileLineNumber;

    return operator OSstream&();
}


Foam::OSstream& Foam::error::operator()
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


Foam::error::operator Foam::OSstream&()
{
    if (!messageStreamPtr_->good())
    {
        Perr<< nl
            << "error::operator OSstream&() : error stream has failed"
            << endl;
        abort();
    }

    return *messageStreamPtr_;
}


Foam::error::operator Foam::dictionary() const
{
    dictionary errDict;

    string oneLineMessage(message());
    oneLineMessage.replaceAll("\n", " ");

    errDict.add("type", word("Foam::error"));
    errDict.add("message", oneLineMessage);
    errDict.add("function", functionName());
    errDict.add("sourceFile", sourceFileName());
    errDict.add("sourceFileLineNumber", sourceFileLineNumber());

    return errDict;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::error::exitOrAbort(const int errNo, const bool isAbort)
{
    if (!throwing_ && JobInfo::constructed)
    {
        jobInfo.add("FatalError", operator dictionary());
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
        error errorException(*this);

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

Foam::string Foam::error::message() const
{
    return messageStreamPtr_->str();
}


void Foam::error::clear() const
{
    return messageStreamPtr_->reset();
}


void Foam::error::exit(const int errNo)
{
    exitOrAbort(errNo, false);
}


void Foam::error::abort()
{
    exitOrAbort(1, true);
}


void Foam::error::write(Ostream& os, const bool includeTitle) const
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


    const label lineNo = sourceFileLineNumber();

    if (error::level >= 2 && lineNo && !functionName().empty())
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

Foam::Ostream& Foam::operator<<(Ostream& os, const error& err)
{
    err.write(os);

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error definitions

Foam::error Foam::FatalError("--> FOAM FATAL ERROR: ");


// ************************************************************************* //
