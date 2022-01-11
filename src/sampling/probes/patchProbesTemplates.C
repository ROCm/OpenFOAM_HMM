/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "patchProbes.H"
#include "volFields.H"
#include "IOmanip.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::patchProbes::sampleAndWrite
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
)
{
    Field<Type> values(sample(vField));

    if (Pstream::master())
    {
        unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& probeStream = *probeFilePtrs_[vField.name()];

        probeStream
            << setw(w)
            << vField.time().timeOutputValue();

        for (const auto& v : values)
        {
            probeStream << ' ' << setw(w) << v;
        }
        probeStream << endl;
    }

    const word& fieldName = vField.name();
    this->setResult("average(" + fieldName + ")", average(values));
    this->setResult("min(" + fieldName + ")", min(values));
    this->setResult("max(" + fieldName + ")", max(values));
    this->setResult("size(" + fieldName + ")", values.size());
}


template<class Type>
void Foam::patchProbes::sampleAndWrite
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
)
{
    Field<Type> values(sample(sField));

    if (Pstream::master())
    {
        unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& probeStream = *probeFilePtrs_[sField.name()];

        probeStream
            << setw(w)
            << sField.time().timeOutputValue();

        for (const auto& v : values)
        {
            probeStream << ' ' << setw(w) << v;
        }
        probeStream << endl;
    }

    const word& fieldName = sField.name();
    this->setResult("average(" + fieldName + ")", average(values));
    this->setResult("min(" + fieldName + ")", min(values));
    this->setResult("max(" + fieldName + ")", max(values));
    this->setResult("size(" + fieldName + ")", values.size());
}


template<class Type>
void Foam::patchProbes::sampleAndWrite
(
    const fieldGroup<Type>& fields
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    for (const auto& fieldName : fields)
    {
        if (loadFromFiles_)
        {
            sampleAndWrite
            (
                VolFieldType
                (
                    IOobject
                    (
                        fieldName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                )
            );
        }
        else
        {
            objectRegistry::const_iterator iter = mesh_.find(fieldName);

            if (iter.found() && iter()->type() == VolFieldType::typeName)
            {
                sampleAndWrite(mesh_.lookupObject<VolFieldType>(fieldName));
            }
        }
    }
}


template<class Type>
void Foam::patchProbes::sampleAndWriteSurfaceFields
(
    const fieldGroup<Type>& fields
)
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    for (const auto& fieldName : fields)
    {
        if (loadFromFiles_)
        {
            sampleAndWrite
            (
                SurfaceFieldType
                (
                    IOobject
                    (
                        fieldName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                )
            );
        }
        else
        {
            objectRegistry::const_iterator iter = mesh_.find(fieldName);

            if (iter.found() && iter()->type() == SurfaceFieldType::typeName)
            {
                sampleAndWrite(mesh_.lookupObject<SurfaceFieldType>(fieldName));
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchProbes::sample
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tValues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tValues.ref();

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(*this, probei)
    {
        label facei = faceList_[probei];

        if (facei >= 0)
        {
            label patchi = patches.whichPatch(facei);
            label localFacei = patches[patchi].whichFace(facei);
            values[probei] = vField.boundaryField()[patchi][localFacei];
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchProbes::sample(const word& fieldName) const
{
    return sample
    (
        mesh_.lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            fieldName
        )
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchProbes::sample
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tValues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tValues.ref();

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(*this, probei)
    {
        label facei = faceList_[probei];

        if (facei >= 0)
        {
            label patchi = patches.whichPatch(facei);
            label localFacei = patches[patchi].whichFace(facei);
            values[probei] = sField.boundaryField()[patchi][localFacei];
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


// ************************************************************************* //
