/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "probes.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOobjectList.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::probes::clearFieldGroups()
{
    scalarFields_.clear();
    vectorFields_.clear();
    sphericalTensorFields_.clear();
    symmTensorFields_.clear();
    tensorFields_.clear();

    surfaceScalarFields_.clear();
    surfaceVectorFields_.clear();
    surfaceSphericalTensorFields_.clear();
    surfaceSymmTensorFields_.clear();
    surfaceTensorFields_.clear();
}


Foam::label Foam::probes::classifyFields()
{
    label nFields = 0;
    clearFieldGroups();

    HashTable<wordHashSet> available =
    (
        loadFromFiles_
      ? IOobjectList(mesh_, mesh_.time().timeName()).classes(fieldSelection_)
      : mesh_.classes(fieldSelection_)
    );

    forAllConstIters(available, iter)
    {
        const word& fieldType = iter.key();
        const wordList fieldNames = iter.val().sortedToc();

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
        else if (fieldType == surfaceScalarField::typeName)
        {
            surfaceScalarFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == surfaceVectorField::typeName)
        {
            surfaceVectorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == surfaceSphericalTensorField::typeName)
        {
            surfaceSphericalTensorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == surfaceSymmTensorField::typeName)
        {
            surfaceSymmTensorFields_.append(fieldNames);
            nFields += count;
        }
        else if (fieldType == surfaceTensorField::typeName)
        {
            surfaceTensorFields_.append(fieldNames);
            nFields += count;
        }
    }

    return nFields;
}


// ************************************************************************* //
