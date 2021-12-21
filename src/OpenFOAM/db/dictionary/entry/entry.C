/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "entry.H"
#include "dictionary.H"
#include "StringStream.H"
#include "JobInfo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::entry::disableFunctionEntries
(
    Foam::debug::infoSwitch("disableFunctionEntries", 0)
);


Foam::entry::inputMode Foam::entry::globalInputMode = inputMode::MERGE;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::entry::reportReadWarning
(
    const IOstream& is,
    const std::string& msg
)
{
    std::cerr
        << "--> FOAM Warning :\n"
        << "    Reading \"" << is.relativeName()
        << "\" at line " << is.lineNumber() << '\n'
        << "    " << msg << std::endl;
}


void Foam::entry::resetInputMode()
{
    globalInputMode = inputMode::MERGE;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entry::entry(const keyType& keyword)
:
    IDLList<entry>::link(),
    keyword_(keyword)
{}


Foam::entry::entry(const entry& e)
:
    IDLList<entry>::link(),
    keyword_(e.keyword_)
{}


Foam::autoPtr<Foam::entry> Foam::entry::clone() const
{
    return clone(dictionary::null);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::entry::raiseBadInput(const ITstream& is) const
{
    const word& keyword = keyword_;

    // Can use FatalIOError instead of SafeFatalIOError
    // since predicate checks are not used at the earliest stages
    FatalIOError
    (
        "",                 // functionName
        "",                 // sourceFileName
        0,                  // sourceFileLineNumber
        this->relativeName(), // ioFileName
        is.lineNumber()     // ioStartLineNumber
    )
        << "Entry '" << keyword << "' with invalid input" << nl << nl
        << exit(FatalIOError);
}


void Foam::entry::checkITstream(const ITstream& is) const
{
    const word& keyword = keyword_;

    if (is.nRemainingTokens())
    {
        const label remaining = is.nRemainingTokens();

        // Similar to SafeFatalIOError
        if (JobInfo::constructed)
        {
            OSstream& err =
                FatalIOError
                (
                    "",                 // functionName
                    "",                 // sourceFileName
                    0,                  // sourceFileLineNumber
                    this->relativeName(), // ioFileName
                    is.lineNumber()     // ioStartLineNumber
                );

            err << "Entry '" << keyword << "' has "
                << remaining << " excess tokens in stream" << nl << nl
                << "    ";
            is.writeList(err, 0);

            err << exit(FatalIOError);
        }
        else
        {
            std::cerr
                << nl
                << "--> FOAM FATAL IO ERROR:" << nl;

            std::cerr
                << "Entry '" << keyword << "' has "
                << remaining << " excess tokens in stream" << nl << nl;

            std::cerr
                << "file: " << this->relativeName()
                << " at line " << is.lineNumber() << '.' << nl
                << std::endl;

            std::exit(1);
        }
    }
    else if (!is.size())
    {
        // Similar to SafeFatalIOError
        if (JobInfo::constructed)
        {
            FatalIOError
            (
                "",                 // functionName
                "",                 // sourceFileName
                0,                  // sourceFileLineNumber
                this->relativeName(), // ioFileName
                is.lineNumber()     // ioStartLineNumber
            )
                << "Entry '" << keyword
                << "' had no tokens in stream" << nl << nl
                << exit(FatalIOError);
        }
        else
        {
            std::cerr
                << nl
                << "--> FOAM FATAL IO ERROR:" << nl
                << "Entry '" << keyword
                << "' had no tokens in stream" << nl << nl;

            std::cerr
                << "file: " << this->relativeName()
                << " at line " << is.lineNumber() << '.' << nl
                << std::endl;

            std::exit(1);
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::entry::operator=(const entry& e)
{
    if (this == &e)
    {
        return;  // Self-assignment is a no-op
    }

    keyword_ = e.keyword_;
}


bool Foam::entry::operator==(const entry& e) const
{
    if (this == &e)
    {
        return true;
    }
    if (keyword_ != e.keyword_)
    {
        return false;
    }

    // Compare contents (as strings)

    OStringStream oss1;
    oss1 << *this;

    OStringStream oss2;
    oss2 << e;

    return oss1.str() == oss2.str();
}


bool Foam::entry::operator!=(const entry& e) const
{
    return !operator==(e);
}


// ************************************************************************* //
