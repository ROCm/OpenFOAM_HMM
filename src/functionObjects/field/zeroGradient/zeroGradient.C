/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "volFields.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(zeroGradient, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        zeroGradient,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::zeroGradient::checkFormatName(const word& str)
{
    if (str.find("@@") == string::npos)
    {
        WarningInFunction
            << "Bad result naming "
            << "(no '@@' token found), deactivating."
            << nl << endl;

        return false;
    }
    else if (str == "@@")
    {
        WarningInFunction
            << "Bad result naming "
            << "(only a '@@' token found), deactivating."
            << nl
            << endl;

        return false;
    }
    else
    {
        return true;
    }
}


void Foam::functionObjects::zeroGradient::uniqWords(wordReList& lst)
{
    boolList retain(lst.size());
    wordHashSet uniq;
    forAll(lst, i)
    {
        const wordRe& select = lst[i];

        retain[i] =
        (
            select.isPattern()
         || uniq.insert(static_cast<const word&>(select))
        );
    }

    inplaceSubset(retain, lst);
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


void Foam::functionObjects::zeroGradient::process()
{
    results_.clear();

    wordHashSet candidates = subsetStrings(selectFields_, mesh_.names());
    DynamicList<word> missing(selectFields_.size());
    DynamicList<word> ignored(selectFields_.size());

    // check exact matches first
    forAll(selectFields_, i)
    {
        const wordRe& select = selectFields_[i];
        if (!select.isPattern())
        {
            const word& fieldName = static_cast<const word&>(select);

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

    forAllConstIter(wordHashSet, candidates, iter)
    {
        process(iter.key());
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::zeroGradient::~zeroGradient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::zeroGradient::read(const dictionary& dict)
{
    dict.lookup("fields") >> selectFields_;
    uniqWords(selectFields_);

    resultName_ = dict.lookupOrDefault<word>("result", type() + "(@@)");
    return checkFormatName(resultName_);
}


bool Foam::functionObjects::zeroGradient::execute()
{
    results_.clear();
    return true;
}


bool Foam::functionObjects::zeroGradient::write()
{
    // Consistent output order
    const wordList outputList = results_.sortedToc();

    forAll(outputList, i)
    {
        const word& fieldName = outputList[i];

        if (foundObject<regIOobject>(fieldName))
        {
            const regIOobject& io = lookupObject<regIOobject>(fieldName);

            Log << type() << " " << name()
                << " write: writing field " << fieldName << endl;

            io.write();
        }
    }

    return true;
}


// ************************************************************************* //
