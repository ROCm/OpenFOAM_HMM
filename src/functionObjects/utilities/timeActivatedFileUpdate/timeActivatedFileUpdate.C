/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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
#include "Time.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(timeActivatedFileUpdate, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        timeActivatedFileUpdate,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::timeActivatedFileUpdate::updateFile()
{
    modified_ = false;

    label i = lastIndex_;
    while
    (
        i < timeVsFile_.size()-1
     && timeVsFile_[i+1].first() < time_.value()+0.5*time_.deltaTValue()
    )
    {
        i++;
    }

    if (i > lastIndex_)
    {
        Log << nl << type() << ": copying file" << nl << timeVsFile_[i].second()
            << nl << "to:" << nl << fileToUpdate_ << nl << endl;

        if (Pstream::master() || time_.distributed())
        {
            // Slaves do not copy if running non-distributed
            fileName destFile(fileToUpdate_ + Foam::name(pid()));
            cp(timeVsFile_[i].second(), destFile);
            mv(destFile, fileToUpdate_);
        }
        lastIndex_ = i;
        modified_ = true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::timeActivatedFileUpdate::timeActivatedFileUpdate
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    fileToUpdate_("unknown-fileToUpdate"),
    timeVsFile_(),
    lastIndex_(-1),
    modified_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::timeActivatedFileUpdate::~timeActivatedFileUpdate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::timeActivatedFileUpdate::read
(
    const dictionary& dict
)
{
    functionObject::read(dict);

    dict.lookup("fileToUpdate") >> fileToUpdate_;
    dict.lookup("timeVsFile") >> timeVsFile_;

    lastIndex_ = -1;
    fileToUpdate_.expand();

    Info<< type() << " " << name() << " output:" << nl
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

        Info<< "    " << timeVsFile_[i].first() << tab
            << timeVsFile_[i].second() << endl;
    }

    // Copy starting files
    updateFile();

    return true;
}


bool Foam::functionObjects::timeActivatedFileUpdate::execute()
{
    updateFile();

    return true;
}


bool Foam::functionObjects::timeActivatedFileUpdate::write()
{
    return true;
}


bool Foam::functionObjects::timeActivatedFileUpdate::filesModified() const
{
    return modified_;
}


// ************************************************************************* //
