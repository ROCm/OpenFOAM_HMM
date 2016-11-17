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
#include "fvMesh.H"
#include "volFields.H"
#include "globalIndex.H"
#include "zeroGradientFvPatchField.H"


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

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            io,
            df.mesh(),
            dimensioned<Type>("0", df.dimensions(), Zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    tvf.ref().primitiveFieldRef() = df;
    tvf.ref().correctBoundaryConditions();

    return tvf;
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
        tmp<GeometricField<Type, fvPatchField, volMesh>> tfld
        (
            subsetter_.interpolate(vf)
        );
        tfld.ref().checkOut();
        tfld.ref().rename(vf.name());
        return tfld;
    }
    else
    {
        return vf;
    }
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
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf =
        zeroGradientField<Type>(df);

    if (subsetter_.hasSubMesh())
    {
        return interpolate<Type>(tvf());
    }
    else
    {
        return tvf;
    }
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
        tmp<GeoField> subFld = subsetter_.interpolate(fld);
        subFld.ref().checkOut();
        subFld.ref().rename(fld.name());
        return subFld;
    }
    else
    {
        return fld;
    }
}


// ************************************************************************* //
