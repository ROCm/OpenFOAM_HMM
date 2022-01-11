/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
#include "IOmanip.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T>
class isNotEqOp
{
public:

    void operator()(T& x, const T& y) const
    {
        const T unsetVal(-VGREAT*pTraits<T>::one);

        if (x != unsetVal)
        {
            // Keep x.

            // Note: should check for y != unsetVal but multiple sample cells
            // already handled in read().
        }
        else
        {
            // x is not set. y might be.
            x = y;
        }
    }
};

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::probes::sampleAndWrite
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
)
{
    Field<Type> values(sample(vField));

    if (Pstream::master())
    {
        unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& os = *probeFilePtrs_[vField.name()];

        os  << setw(w) << vField.time().timeOutputValue();

        forAll(values, probei)
        {
            if (includeOutOfBounds_ || processor_[probei] != -1)
            {
                os  << ' ' << setw(w) << values[probei];
            }
        }
        os  << endl;
    }

    const word& fieldName = vField.name();
    this->setResult("average(" + fieldName + ")", average(values));
    this->setResult("min(" + fieldName + ")", min(values));
    this->setResult("max(" + fieldName + ")", max(values));
    this->setResult("size(" + fieldName + ")", values.size());
}


template<class Type>
void Foam::probes::sampleAndWrite
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
)
{
    Field<Type> values(sample(sField));

    if (Pstream::master())
    {
        unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& os = *probeFilePtrs_[sField.name()];

        os  << setw(w) << sField.time().timeOutputValue();

        forAll(values, probei)
        {
            if (includeOutOfBounds_ || processor_[probei] != -1)
            {
                os  << ' ' << setw(w) << values[probei];
            }
        }
        os  << endl;
    }

    const word& fieldName = sField.name();
    this->setResult("average(" + fieldName + ")", average(values));
    this->setResult("min(" + fieldName + ")", min(values));
    this->setResult("max(" + fieldName + ")", max(values));
    this->setResult("size(" + fieldName + ")", values.size());
}


template<class Type>
void Foam::probes::sampleAndWrite(const fieldGroup<Type>& fields)
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
void Foam::probes::sampleAndWriteSurfaceFields(const fieldGroup<Type>& fields)
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
Foam::probes::sample
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tValues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tValues.ref();

    if (fixedLocations_)
    {
        autoPtr<interpolation<Type>> interpolator
        (
            interpolation<Type>::New(interpolationScheme_, vField)
        );

        forAll(*this, probei)
        {
            if (elementList_[probei] >= 0)
            {
                const vector& position = operator[](probei);

                values[probei] = interpolator().interpolate
                (
                    position,
                    elementList_[probei],
                    -1
                );
            }
        }
    }
    else
    {
        forAll(*this, probei)
        {
            if (elementList_[probei] >= 0)
            {
                values[probei] = vField[elementList_[probei]];
            }
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sample(const word& fieldName) const
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
Foam::probes::sample
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tValues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tValues.ref();

    forAll(*this, probei)
    {
        if (faceList_[probei] >= 0)
        {
            values[probei] = sField[faceList_[probei]];
        }
    }

    Pstream::listCombineGather(values, isNotEqOp<Type>());
    Pstream::listCombineScatter(values);

    return tValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sampleSurfaceFields(const word& fieldName) const
{
    return sample
    (
        mesh_.lookupObject<GeometricField<Type, fvsPatchField, surfaceMesh>>
        (
            fieldName
        )
    );
}

// ************************************************************************* //
