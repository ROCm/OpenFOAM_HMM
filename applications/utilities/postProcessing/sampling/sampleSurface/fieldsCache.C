/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fieldsCache.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fieldsCache<Type>::fieldsCache()
:
    HashPtrTable<GeometricField<Type, fvPatchField, volMesh> >(),
    pointFields_(),
    interpolators_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>&
Foam::fieldsCache<Type>::pointField
(
    const word& name,
    const volPointInterpolation& pInterp
) const
{
    if (!pointFields_.found(name))
    {
        const GeometricField<Type, fvPatchField, volMesh>& vField =
            *this->operator[](name);

        tmp<GeometricField<Type, pointPatchField, pointMesh> > tptField =
            pInterp.interpolate(vField);

        GeometricField<Type, pointPatchField, pointMesh>* ptFieldPtr =
            tptField.ptr();

        pointFields_.insert(name, ptFieldPtr);

        return *ptFieldPtr;
    }
    else
    {
        return *pointFields_[name];
    }
}


template<class Type>
const Foam::interpolation<Type>& Foam::fieldsCache<Type>::interpolator
(
    const word& name,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    if (!interpolators_.found(name))
    {
        const GeometricField<Type, fvPatchField, volMesh>& vField =
            *this->operator[](name);

        interpolation<Type>* interpolatorPtr =
            interpolation<Type>::New
            (
                interpolationSchemes,
                pInterp,
                vField
            ).ptr();

        interpolators_.insert(name, interpolatorPtr);

        return *interpolatorPtr;
    }
    else
    {
        return *interpolators_[name];
    }
}


// ************************************************************************* //
