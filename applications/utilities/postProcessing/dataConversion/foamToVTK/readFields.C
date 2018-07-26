/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "readFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class GeoField>
label readFields
(
    const fvMeshSubsetProxy& proxy,
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const wordHashSet& selectedFields,
    PtrList<const GeoField>& fields
)
{
    // Available fields of type GeoField, sorted order
    const wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.sortedNames(GeoField::typeName)
      : objects.sortedNames(GeoField::typeName, selectedFields)
    );

    // Construct the fields
    fields.resize(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        fields.set
        (
            nFields++,
            proxy.interpolate
            (
                GeoField(*(objects[fieldName]), mesh)
            ).ptr()
        );
    }

    return nFields;
}


template<class GeoField>
label readFields
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const wordHashSet& selectedFields,
    PtrList<const GeoField>& fields
)
{
    // Available fields of type GeoField, sorted order
    const wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.sortedNames(GeoField::typeName)
      : objects.sortedNames(GeoField::typeName, selectedFields)
    );

    // Construct the fields
    fields.resize(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        fields.set
        (
            nFields++,
            new GeoField(*(objects[fieldName]), mesh)
        );
    }

    return nFields;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
