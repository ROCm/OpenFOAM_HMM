/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "zeroGradient.H"
#include "stringListOps.H"
#include "volFields.H"
#include "dictionary.H"
#include "wordRes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(zeroGradient, 0);
    addToRunTimeSelectionTable(functionObject, zeroGradient, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::zeroGradient::checkFormatName
(
    const std::string& str
)
{
    if (std::string::npos == str.find("@@"))
    {
        WarningInFunction
            << "Bad result naming (no '@@' token found)."
            << nl << endl;

        return false;
    }
    else if (str == "@@")
    {
        WarningInFunction
            << "Bad result naming (only a '@@' token found)."
            << nl << endl;

        return false;
    }

    return true;
}


int Foam::functionObjects::zeroGradient::process(const word& fieldName)
{
    int state = 0;
    apply<scalar>(fieldName, state);
    apply<vector>(fieldName, state);
    apply<sphericalTensor>(fieldName, state);
    apply<symmTensor>(fieldName, state);
    apply<tensor>(fieldName, state);

    return state;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::zeroGradient::zeroGradient
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    selectFields_(),
    resultName_(string::null),
    results_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::zeroGradient::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readEntry("fields", selectFields_);
    selectFields_.uniq();

    Info<< type() << " fields: " << selectFields_ << nl;

    resultName_ =
        dict.getOrDefault<word>("result", scopedName(type() + "(@@)"));

    // Require '@@' token for result, unless a single (non-regex) source field
    return
    (
        (selectFields_.size() == 1 && selectFields_.first().isLiteral())
     || checkFormatName(resultName_)
    );
}


bool Foam::functionObjects::zeroGradient::execute()
{
    results_.clear();

    wordHashSet candidates(subsetStrings(selectFields_, mesh_.names()));
    DynamicList<word> missing(selectFields_.size());
    DynamicList<word> ignored(selectFields_.size());

    // Check exact matches first
    for (const wordRe& select : selectFields_)
    {
        if (select.isLiteral())
        {
            const word& fieldName = select;

            if (!candidates.erase(fieldName))
            {
                missing.append(fieldName);
            }
            else if (process(fieldName) < 1)
            {
                ignored.append(fieldName);
            }
        }
    }

    for (const word& fieldName : candidates)
    {
        process(fieldName);
    }

    if (missing.size())
    {
        WarningInFunction
            << "Missing field " << missing << endl;
    }
    if (ignored.size())
    {
        WarningInFunction
            << "Unprocessed field " << ignored << endl;
    }

    return true;
}


bool Foam::functionObjects::zeroGradient::write()
{
    if (results_.size())
    {
        Log << type() << ' ' << name() << " write:" << endl;
    }

    // Consistent output order
    for (const word& fieldName : results_.sortedToc())
    {
        const regIOobject* ioptr = findObject<regIOobject>(fieldName);

        if (ioptr)
        {
            Log << "    " << fieldName << endl;

            ioptr->write();
        }
    }

    return true;
}


// ************************************************************************* //
