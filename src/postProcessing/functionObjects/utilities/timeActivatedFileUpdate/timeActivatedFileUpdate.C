/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "timeActivatedFileUpdate.H"
#include "objectRegistry.H"
#include "Time.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeActivatedFileUpdate, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::timeActivatedFileUpdate::updateFile()
{
    label i = lastIndex_;
    while
    (
        i < timeVsFile_.size()-1
     && timeVsFile_[i+1].first() < obr_.time().value()
    )
    {
        i++;
    }

    if (i > lastIndex_)
    {
        if (log_) Info
            << nl << type() << ": copying file" << nl << timeVsFile_[i].second()
            << nl << "to:" << nl << fileToUpdate_ << nl << endl;

        cp(timeVsFile_[i].second(), fileToUpdate_);
        lastIndex_ = i;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeActivatedFileUpdate::timeActivatedFileUpdate
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    fileToUpdate_("unknown-fileToUpdate"),
    timeVsFile_(),
    lastIndex_(-1),
    log_(true)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeActivatedFileUpdate::~timeActivatedFileUpdate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeActivatedFileUpdate::read(const dictionary& dict)
{
    if (active_)
    {
        log_.readIfPresent("log", dict);

        dict.lookup("fileToUpdate") >> fileToUpdate_;
        dict.lookup("timeVsFile") >> timeVsFile_;

        lastIndex_ = -1;
        fileToUpdate_.expand();

        if (log_) Info
            << type() << " " << name_ << " output:" << nl
            << "    time vs file list:" << endl;

        forAll(timeVsFile_, i)
        {
            timeVsFile_[i].second() = timeVsFile_[i].second().expand();
            if (!isFile(timeVsFile_[i].second()))
            {
                FatalErrorInFunction
                    << "File: " << timeVsFile_[i].second() << " not found"
                    << nl << exit(FatalError);
            }

            if (log_) Info
                << "    " << timeVsFile_[i].first() << tab
                << timeVsFile_[i].second() << endl;
        }
        if (log_) Info<< endl;

        updateFile();
    }
}


void Foam::timeActivatedFileUpdate::execute()
{
    if (active_)
    {
        updateFile();
    }
}


void Foam::timeActivatedFileUpdate::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::timeActivatedFileUpdate::timeSet()
{
    // Do nothing
}


void Foam::timeActivatedFileUpdate::write()
{
    // Do nothing
}


// ************************************************************************* //
