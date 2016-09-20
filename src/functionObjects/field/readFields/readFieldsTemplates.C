/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::readFields::loadField(const word& fieldName) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (obr_.foundObject<VolFieldType>(fieldName))
    {
        DebugInfo
            << "readFields : " << VolFieldType::typeName
            << " " << fieldName << " already in database"
            << endl;
    }
    else if (obr_.foundObject<SurfaceFieldType>(fieldName))
    {
        DebugInfo<< "readFields: " << SurfaceFieldType::typeName
            << " " << fieldName << " already exists in database"
            << " already in database" << endl;
    }
    else
    {
        IOobject fieldHeader
        (
            fieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if
        (
            fieldHeader.typeHeaderOk<vfType>(false)
         && fieldHeader.headerClassName() == VolFieldType::typeName
        )
        {
            // Store field on mesh database
            Log << "    Reading " << fieldName << endl;
            tmp<VolFieldType> tvf(new VolFieldType(fieldHeader, mesh));
            store(tvf, fieldName);
            return true;
        }
        else if
        (
            fieldHeader.typeHeaderOk<sfType>(false)
         && fieldHeader.headerClassName() == SurfaceFieldType::typeName
        )
        {
            // Store field on mesh database
            Log << "    Reading " << fieldName << endl;
            tmp<SurfaceFieldType> tsf(new SurfaceFieldType(fieldHeader, mesh));
            store(tsf, fieldName);
            return true;
        }
    }

    return false;
}


// ************************************************************************* //
