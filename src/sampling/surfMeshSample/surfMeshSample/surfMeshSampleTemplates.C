/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "surfMeshSample.H"
#include "dimensionedType.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::surfMeshSample::sampleOnFaces
(
    const interpolation<Type>& sampler,
    const labelUList& elements,
    const faceList& fcs,
    const pointField& pts
)
{
    const label len = elements.size();

    if (len != fcs.size())
    {
        FatalErrorInFunction
            << "size mismatch: "
            << "sampled elements (" << len
            << ") != faces (" << fcs.size() << ')'
            << exit(FatalError);
    }

    auto tvalues = tmp<Field<Type>>::New(len);
    auto& values = tvalues.ref();

    for (label i=0; i < len; ++i)
    {
        const label celli = elements[i];
        const point pt = fcs[i].centre(pts);

        values[i] = sampler.interpolate(pt, celli);
    }

    return tvalues;
}


template<class Type>
Foam::DimensionedField<Type, Foam::surfGeoMesh>&
Foam::surfMeshSample::getOrCreateSurfField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    const surfMesh& surf  = surface();
    const word& fieldName = vField.name();

    SurfFieldType* ptr = surf.getObjectPtr<SurfFieldType>(fieldName);
    if (!ptr)
    {
        ptr = new SurfFieldType
        (
            IOobject
            (
                fieldName,
                surf.time().timeName(),
                surf,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            surf,
            dimensioned<Type>(vField.dimensions(), Zero)
        );
        ptr->writeOpt() = IOobject::NO_WRITE;

        surf.store(ptr);
    }

    return *ptr;
}


// // Older code for transferring an IOField to a surfField between
// // different registries
// template<class Type>
// bool Foam::surfMeshSample::transferField
// (
//     const objectRegistry& store,
//     const word& fieldName
// ) const
// {
//     typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
//     typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;
//     typedef IOField<Type> TmpFieldType;
//
//     // foundObject includes a type check
//     bool ok =
//     (
//         mesh_.foundObject<VolFieldType>(fieldName)
//      && store.foundObject<TmpFieldType>(fieldName)
//     );
//
//     if (ok)
//     {
//         SurfFieldType& sfield = getOrCreateSurfField
//         (
//             mesh_.lookupObject<VolFieldType>(fieldName)
//         );
//
//         TmpFieldType& iofield =
//             store.lookupObjectRef<TmpFieldType>(fieldName);
//
//         sfield.transfer(iofield);
//         store.checkOut(iofield);
//     }
//
//     return ok;
// }


template<class Type>
Foam::label Foam::surfMeshSample::writeFields
(
    const wordRes& select
) const
{
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;
    const surfMesh& s = surface();

    const wordList fieldNames(s.sortedNames<SurfFieldType>(select));

    for (const word& fieldName : fieldNames)
    {
        s.lookupObject<SurfFieldType>(fieldName).write();
    }

    return fieldNames.size();
}


// ************************************************************************* //
