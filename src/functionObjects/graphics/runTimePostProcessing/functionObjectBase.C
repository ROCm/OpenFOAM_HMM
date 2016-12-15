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

#include "functionObjectBase.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{
    defineTypeNameAndDebug(functionObjectBase, 0);
}
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostPro::functionObjectBase::removeFile
(
    const word& keyword,
    const word& fieldName
)
{
    dictionary dict;
    state_.getObjectDict(functionObjectName_, fieldName, dict);

    fileName fName;
    if (dict.readIfPresent(keyword, fName))
    {
        Foam::rm(fName);
        return true;
    }

    return false;
}


Foam::fileName
Foam::functionObjects::runTimePostPro::functionObjectBase::getFileName
(
    const word& keyword,
    const word& fieldName
) const
{
    dictionary dict;
    state_.getObjectDict(functionObjectName_, fieldName, dict);

    fileName fName(dict.lookupOrDefault(keyword, fileName::null));

    return fName;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectBase::functionObjectBase
(
    const stateFunctionObject& state,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>, word>& colours
)
:
    fieldVisualisationBase(dict, colours),
    state_(state),
    functionObjectName_(""),
    clearObjects_(dict.lookupOrDefault<bool>("clearObjects", false))
{
    dict.lookup("functionObject") >> functionObjectName_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectBase::~functionObjectBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostPro::functionObjectBase::clear()
{
    return clearObjects_;
}


// ************************************************************************* //
