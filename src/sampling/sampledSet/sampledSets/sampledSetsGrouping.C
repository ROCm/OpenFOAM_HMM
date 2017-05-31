/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "sampledSets.H"
#include "volFields.H"
#include "IOobjectList.H"
#include "stringListOps.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::clearFieldGroups()
{
    scalarFields_.clear();
    vectorFields_.clear();
    sphericalTensorFields_.clear();
    symmTensorFields_.clear();
    tensorFields_.clear();
}


Foam::label Foam::sampledSets::classifyFields()
{
    label nFields = 0;
    clearFieldGroups();

    wordList allFields;    // Just needed for warnings
    HashTable<wordHashSet> available;

    if (loadFromFiles_)
    {
        // Check files for a particular time
        IOobjectList objects(mesh_, mesh_.time().timeName());

        allFields = objects.names();
        available = objects.classes(fieldSelection_);
    }
    else
    {
        // Check currently available fields
        allFields = mesh_.names();
        available = mesh_.classes(fieldSelection_);
    }

    DynamicList<label> missed(fieldSelection_.size());

    // Detect missing fields
    forAll(fieldSelection_, i)
    {
        if (findStrings(fieldSelection_[i], allFields).empty())
        {
            missed.append(i);
        }
    }

    if (missed.size())
    {
        WarningInFunction
            << nl
            << "Cannot find "
            << (loadFromFiles_ ? "field file" : "registered field")
            << " matching "
            << UIndirectList<wordRe>(fieldSelection_, missed) << endl;
    }

    forAllConstIters(available, iter)
    {
        const word& fieldType = iter.key();
        const wordList fieldNames = iter.object().sortedToc();

        const label count = fieldNames.size(); // pre-filtered, so non-empty

        if (fieldType == volScalarField::typeName)
        {
            scalarFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == volVectorField::typeName)
        {
            vectorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == volSphericalTensorField::typeName)
        {
            sphericalTensorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == volSymmTensorField::typeName)
        {
            symmTensorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == volTensorField::typeName)
        {
            tensorFields_.append(fieldNames);
            nFields += count;
        }
    }

    return nFields;
}


// ************************************************************************* //
