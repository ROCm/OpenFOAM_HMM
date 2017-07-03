/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "abort.H"
#include "dictionary.H"
#include "error.H"
#include "Time.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(abort, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        abort,
        dictionary
    );
}
}


const Foam::Enum
<
    Foam::Time::stopAtControls
>
Foam::functionObjects::abort::actionNames_
{
    { Time::stopAtControls::saNoWriteNow, "noWriteNow" },
    { Time::stopAtControls::saWriteNow, "writeNow" },
    { Time::stopAtControls::saNextWrite, "nextWrite" },
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::abort::removeFile() const
{
    bool hasAbort = isFile(abortFile_);
    reduce(hasAbort, orOp<bool>());

    if (hasAbort && Pstream::master())
    {
        // Cleanup ABORT file (on master only)
        rm(abortFile_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::abort::abort
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    abortFile_("$FOAM_CASE/" + name),
    action_(Time::stopAtControls::saNextWrite)
{
    abortFile_.expand();
    read(dict);

    // Remove any old files from previous runs
    removeFile();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::abort::~abort()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::abort::read(const dictionary& dict)
{
    functionObject::read(dict);

    action_ = actionNames_.lookupOrDefault
    (
        "action",
        dict,
        Time::stopAtControls::saNextWrite
    );

    if (dict.readIfPresent("file", abortFile_))
    {
        abortFile_.expand();
    }

    return true;
}


bool Foam::functionObjects::abort::execute()
{
    bool hasAbort = isFile(abortFile_);
    reduce(hasAbort, orOp<bool>());

    if (hasAbort)
    {
        switch (action_)
        {
            case Time::saNoWriteNow :
            {
                if (time_.stopAt(action_))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << time_.timeIndex()
                        << "): stop without writing data"
                        << endl;
                }
                break;
            }

            case Time::saWriteNow :
            {
                if (time_.stopAt(action_))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << time_.timeIndex()
                        << "): stop+write data"
                        << endl;
                }
                break;
            }

            case Time::saNextWrite :
            {
                if (time_.stopAt(action_))
                {
                    Info<< "USER REQUESTED ABORT (timeIndex="
                        << time_.timeIndex()
                        << "): stop after next data write"
                        << endl;
                }
                break;
            }

            default:
            {
                // Invalid choices already filtered out by Enum
            }
        }
    }

    return true;
}


bool Foam::functionObjects::abort::write()
{
    return true;
}


bool Foam::functionObjects::abort::end()
{
    removeFile();
    return true;
}


// ************************************************************************* //
