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

#include "meshToMesh.H"

template<class Type>
bool Foam::mapFieldsFO::writeFieldType() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    const fvMesh& mapRegion = mapRegionPtr_();

    wordList fieldNames(obr_.names(FieldType::typeName));
    labelList selected = findStrings(fieldNames_, fieldNames);
    forAll(selected, i)
    {
        const word& fieldName = fieldNames[selected[i]];
        const FieldType& field = obr_.lookupObject<FieldType>(fieldName);

        if (log_) Info << "    " << fieldName;

        IOobject mapRegionIO
        (
            fieldName,
            field.time().timeName(),
            mapRegion
        );

        tmp<FieldType> tfieldMapRegion(interpPtr_->mapTgtToSrc(field));

        if (log_) Info<< ": interpolated";

        FieldType fieldMapRegion(mapRegionIO, tfieldMapRegion);
        fieldMapRegion.write();

        if (log_) Info<< " and written" << nl;
    }

    return selected.size() > 0;
}


// ************************************************************************* //
