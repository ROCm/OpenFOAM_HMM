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

#include "meshSubsetHelper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::meshSubsetHelper::zeroGradientField
(
    const typename GeometricField
    <
        Type,
        fvPatchField,
        volMesh
    >::Internal& df
)
{
    IOobject io(df);
    io.readOpt()  = IOobject::NO_READ;
    io.writeOpt() = IOobject::NO_WRITE;
    io.registerObject() = false;

    auto tfield = tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        io,
        df.mesh(),
        dimensioned<Type>(df.dimensions(), Zero),
        zeroGradientFvPatchField<Type>::typeName
    );
    tfield.ref().primitiveFieldRef() = df;
    tfield.ref().oriented() = df.oriented();
    tfield.ref().correctBoundaryConditions();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::meshSubsetHelper::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (subsetter_.hasSubMesh())
    {
        auto tfield(subsetter_.interpolate(vf));

        tfield.ref().checkOut();
        tfield.ref().rename(vf.name());
        return tfield;
    }

    return vf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::meshSubsetHelper::interpolate
(
    const typename GeometricField
    <
        Type,
        fvPatchField,
        volMesh
    >::Internal& df
) const
{
    auto tfield = zeroGradientField<Type>(df);

    if (subsetter_.hasSubMesh())
    {
        return interpolate<Type>(tfield());
    }

    return tfield;
}


template<class GeoField>
Foam::tmp<GeoField>
Foam::meshSubsetHelper::interpolate
(
    const GeoField& fld
) const
{
    if (subsetter_.hasSubMesh())
    {
        tmp<GeoField> tfield = subsetter_.interpolate(fld);
        tfield.ref().checkOut();
        tfield.ref().rename(fld.name());
        return tfield;
    }

    return fld;
}


// ************************************************************************* //
