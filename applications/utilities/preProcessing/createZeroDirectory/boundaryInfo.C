/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "boundaryInfo.H"
#include "Time.H"
#include "IOPtrList.H"
#include "polyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::boundaryInfo::boundaryInfo(const Time& runTime, const word& regionName)
:
    names_(),
    types_(),
    constraint_(),
    groups_(),
    allGroupNames_()
{
    // read the mesh boundaries for the current region
    Info<< "    Reading mesh boundaries" << endl;

    const_cast<word&>(IOPtrList<entry>::typeName) = word::null;
    IOPtrList<entry> boundaryPatchList
    (
        IOobject
        (
            "boundary",
            runTime.constant(),
            regionName/polyMesh::meshSubDir,
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    names_.setSize(boundaryPatchList.size());
    types_.setSize(boundaryPatchList.size());
    constraint_.setSize(boundaryPatchList.size(), false);
    groups_.setSize(boundaryPatchList.size());


    forAll(boundaryPatchList, patchI)
    {
        const dictionary& dict = boundaryPatchList[patchI].dict();

        names_[patchI] = dict.dictName();
        dict.lookup("type") >> types_[patchI];
        if (polyPatch::constraintType(types_[patchI]))
        {
            constraint_[patchI] = true;
        }

        if (dict.found("inGroups"))
        {
            dict.lookup("inGroups") >> groups_[patchI];
            allGroupNames_.insert(groups_[patchI]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::boundaryInfo::names() const
{
    return names_;
}


const Foam::wordList& Foam::boundaryInfo::types() const
{
    return types_;
}


const Foam::boolList& Foam::boundaryInfo::constraint() const
{
    return constraint_;
}


const Foam::List<Foam::wordList>& Foam::boundaryInfo::groups() const
{
    return groups_;
}


const Foam::wordHashSet& Foam::boundaryInfo::allGroupNames() const
{
    return allGroupNames_;
}


// ************************************************************************* //
