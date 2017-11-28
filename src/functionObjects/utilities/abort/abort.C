/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// file-scope
// Long description for the action name
namespace Foam
{
static std::string longDescription(const Time::stopAtControls ctrl)
{
    switch (ctrl)
    {
        case Foam::Time::saNoWriteNow :
        {
            return "stop without writing data";
            break;
        }

        case Time::saWriteNow :
        {
            return "stop and write data";
            break;
        }

        case Time::saNextWrite :
        {
            return "stop after next data write";
            break;
        }

        default:
        {
            // Invalid choices already filtered out by Enum
            return "abort";
            break;
        }
    }
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
    action_(Time::stopAtControls::saNextWrite),
    triggered_(false)
{
    abortFile_.expand();
    read(dict);

    // Cleanup old files from previous runs
    if (Pstream::master())
    {
        Foam::rm(abortFile_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::abort::read(const dictionary& dict)
{
    functionObject::read(dict);

    if (dict.readIfPresent("file", abortFile_))
    {
        abortFile_.expand();
    }

    const auto oldAction = action_;

    action_ = Time::stopAtControlNames.lookupOrDefault
    (
        "action",
        dict,
        Time::stopAtControls::saNextWrite
    );

    // User can change action and re-trigger the abort.
    // eg, they had nextWrite, but actually wanted writeNow.
    if (oldAction != action_)
    {
        triggered_ = false;
    }

    Info<< type() << " activated ("
        << longDescription(action_).c_str() <<")" << nl
        << "    File: " << abortFile_ << endl;

    return true;
}


bool Foam::functionObjects::abort::execute()
{
    // If it has been triggered (eg, nextWrite) don't need to check it again
    if (!triggered_)
    {
        bool hasAbort = (Pstream::master() && isFile(abortFile_));
        Pstream::scatter(hasAbort);

        if (hasAbort)
        {
            triggered_ = time_.stopAt(action_);

            if (triggered_)
            {
                Info<< "USER REQUESTED ABORT (timeIndex="
                    << time_.timeIndex()
                    << "): " << longDescription(action_).c_str()
                    << endl;
            }

            Pstream::scatter(triggered_);
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
    // Cleanup ABORT file
    if (Pstream::master())
    {
        Foam::rm(abortFile_);
    }

    return true;
}


// ************************************************************************* //
