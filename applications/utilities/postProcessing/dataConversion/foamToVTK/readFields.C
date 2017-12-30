/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "readFields.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class GeoField>
label readFields
(
    const meshSubsetHelper& helper,
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    PtrList<const GeoField>& fields
)
{
    label nFields = 0;

    // Available fields of type GeomField
    const wordList fieldNames = objects.sortedNames(GeoField::typeName);

    fields.setSize(fieldNames.size());

    // Construct the fields
    for (const word& fieldName : fieldNames)
    {
        if (selectedFields.empty() || selectedFields.found(fieldName))
        {
            fields.set
            (
                nFields++,
                helper.interpolate
                (
                    GeoField(*(objects[fieldName]), mesh)
                ).ptr()
            );
        }
    }

    fields.setSize(nFields);

    return nFields;
}


template<class GeoField>
void readFields
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    PtrList<const GeoField>& fields
)
{
    // Search list of objects for fields of type GeomField
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

    // Construct the fields
    fields.setSize(fieldObjects.size());
    label nFields = 0;

    forAllIters(fieldObjects, iter)
    {
        if (selectedFields.empty() || selectedFields.found(iter()->name()))
        {
            fields.set
            (
                nFields++,
                new GeoField(*iter(), mesh)
            );
        }
    }

    fields.setSize(nFields);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
