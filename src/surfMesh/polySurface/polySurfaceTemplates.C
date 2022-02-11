/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "polySurface.H"
#include "polySurfaceFields.H"
#include "polySurfacePointFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoMeshType>
const Foam::regIOobject* Foam::polySurface::findFieldObject
(
    const word& fieldName
) const
{
    // Face Data first (main registry)

    const objectRegistry& obr = *this;

    const auto* ioptr = obr.cfindObject<regIOobject>(fieldName);

    if (ioptr)
    {
        return ioptr;
    }

    forAllConstIters(obr, iter)
    {
        const objectRegistry* subreg = isA<objectRegistry>(iter.val());

        if (subreg && (ioptr = subreg->cfindObject<regIOobject>(fieldName)))
        {
            return ioptr;
        }
    }

    return ioptr;
}


template<class GeoMeshType>
const Foam::objectRegistry* Foam::polySurface::whichRegistry
(
    const word& fieldName
) const
{
    // Face Data first (main registry)

    const objectRegistry& obr = *this;

    if (obr.found(fieldName))
    {
        return this;
    }

    forAllConstIters(obr, iter)
    {
        const objectRegistry* subreg = isA<objectRegistry>(iter.val());

        if (subreg && subreg->found(fieldName))
        {
            return subreg;
        }
    }

    return nullptr;
}


template<class Type, class GeoMeshType>
void Foam::polySurface::storeField
(
    const word& fieldName,
    const dimensionSet& dims,
    const Field<Type>& values
)
{
    // Force creates field database if needed.
    const objectRegistry& fieldDb = this->fieldData<GeoMeshType>();

    auto* dimfield =
        fieldDb.getObjectPtr<DimensionedField<Type, GeoMeshType>>(fieldName);

    if (dimfield)
    {
        dimfield->dimensions().reset(dims);  // Dimensions may have changed
        dimfield->field() = values;
    }
    else
    {
        dimfield = new DimensionedField<Type, GeoMeshType>
        (
            IOobject
            (
                fieldName,
                fieldDb,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this,
            dims,
            values
        );

        dimfield->store();
    }
}


template<class Type, class GeoMeshType>
void Foam::polySurface::storeField
(
    const word& fieldName,
    const dimensionSet& dims,
    Field<Type>&& values
)
{
    // Force creates field database if needed.
    const objectRegistry& fieldDb = this->fieldData<GeoMeshType>();

    auto* dimfield =
        fieldDb.getObjectPtr<DimensionedField<Type, GeoMeshType>>(fieldName);

    if (dimfield)
    {
        dimfield->dimensions().reset(dims);  // Dimensions may have changed
        dimfield->field() = std::move(values);
    }
    else
    {
        dimfield = new DimensionedField<Type, GeoMeshType>
        (
            IOobject
            (
                fieldName,
                fieldDb,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this,
            dims,
            std::move(values)
        );

        dimfield->store();
    }
}


// ************************************************************************* //
