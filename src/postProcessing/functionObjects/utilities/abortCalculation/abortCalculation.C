/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "abortCalculation.H"
#include "dictionary.H"
#include "error.H"
#include "Time.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(abortCalculation, 0);
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum<Foam::abortCalculation::actionType, 3>::names[] =
{
    "noWriteNow",
    "writeNow",
    "nextWrite"
};

const Foam::NamedEnum<Foam::abortCalculation::actionType, 3>
    Foam::abortCalculation::actionTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::abortCalculation::removeFile() const
{
    if (isFile(abortFile_))
    {
        rm(abortFile_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::abortCalculation::abortCalculation
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    abortFile_("$FOAM_CASE/" + name),
    action_(nextWrite)
{
    abortFile_.expand();
    read(dict);

    // remove any old files from previous runs
    removeFile();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::abortCalculation::~abortCalculation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::abortCalculation::read(const dictionary& dict)
{
    word actionName;

    if (dict.found("action"))
    {
        action_ = actionTypeNames_.read
        (
            dict.lookup("action")
        );
    }
    else
    {
        action_ = nextWrite;
    }

    if (dict.readIfPresent("fileName", abortFile_))
    {
        abortFile_.expand();
    }
}


void Foam::abortCalculation::execute()
{
    if (isFile(abortFile_))
    {
        switch (action_)
        {
            case noWriteNow :
                obr_.time().stopAt(Time::saNoWriteNow);
                Info<< "user requested abort - "
                       "stop immediately without writing data" << endl;
                break;

            case writeNow :
                obr_.time().stopAt(Time::saWriteNow);
                Info<< "user requested abort - "
                       "stop immediately with writing data" << endl;
                break;

            case nextWrite :
                obr_.time().stopAt(Time::saNextWrite);
                Info<< "user requested abort - "
                       "stop after next data write" << endl;
                break;
        }
    }
}


void Foam::abortCalculation::end()
{
    removeFile();
}


void Foam::abortCalculation::write()
{
    execute();
}


// ************************************************************************* //
