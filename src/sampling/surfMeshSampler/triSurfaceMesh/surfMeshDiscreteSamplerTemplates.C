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

#include "surfMeshDiscreteSampler.H"
#include "dimensionedType.H"
#include "error.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::surfMeshDiscreteSampler::sampleType
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const polyMesh& mesh = SurfaceSource::mesh();

    if (!mesh.foundObject<VolFieldType>(fieldName))
    {
        return false;
    }

    const VolFieldType& fld = mesh.lookupObject<VolFieldType>(fieldName);

    getOrCreateSurfField<Type>(fld).field()
        = SurfaceSource::sampleField(fld);

    return true;
}


// template<class Type>
// Foam::tmp<Foam::DimensionedField<Type, Foam::surfGeoMesh>>
// Foam::surfMeshDiscreteSampler::sampleField
// (
//     const GeometricField<Type, fvPatchField, volMesh>& fld
// ) const
// {
//     typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;
//
//     tmp<Field<Type>> tfield = SurfaceSource::sampleField(fld);
//     SurfFieldType& result = getOrCreateSurfField<Type>(fld);
//
//     // need to verify if this will be needed (in the future)
//     const surfMesh& s = surface();
//     if (result.size() != s.size())
//     {
//         // maybe resampling changed the surfMesh,
//         // but the existing surfField wasn't updated
//
//         result.setSize(s.size(), Zero);
//     }
//
//     if (result.size() != sampleElements().size())
//     {
//         FatalErrorInFunction
//             << "mismatch in field/mesh sizes "
//             << result.name() << nl
//             << " field has " << result.size()
//             << " but sampleElements has " << sampleElements().size() << nl
//             << endl
//             << exit(FatalError);
//
//         Info<< "WARNING: "
//             << endl;
//     }
//
//     result = tfield;
//     return result;
// }


// ************************************************************************* //
