/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "GeometricSymmTensorField.H"
#include "symmTensorFieldField.H"

#define TEMPLATE template<template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.C"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Cmpt, template<class> class PatchField, class GeoMesh>
void Foam::zip
(
    GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& result,
    const GeometricField<Cmpt, PatchField, GeoMesh>& xx,
    const GeometricField<Cmpt, PatchField, GeoMesh>& xy,
    const GeometricField<Cmpt, PatchField, GeoMesh>& xz,
    const GeometricField<Cmpt, PatchField, GeoMesh>& yy,
    const GeometricField<Cmpt, PatchField, GeoMesh>& yz,
    const GeometricField<Cmpt, PatchField, GeoMesh>& zz
)
{
    Foam::zip
    (
        result.primitiveFieldRef(),
        xx.primitiveField(), xy.primitiveField(), xz.primitiveField(),
        yy.primitiveField(), yz.primitiveField(),
        zz.primitiveField()
    );

    Foam::zip
    (
        result.boundaryFieldRef(),
        xx.boundaryField(), xy.boundaryField(), xz.boundaryField(),
        yy.boundaryField(), yz.boundaryField(),
        zz.boundaryField()
    );
}


template<class Cmpt, template<class> class PatchField, class GeoMesh>
void Foam::unzip
(
    const GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Cmpt, PatchField, GeoMesh>& xx,
    GeometricField<Cmpt, PatchField, GeoMesh>& xy,
    GeometricField<Cmpt, PatchField, GeoMesh>& xz,
    GeometricField<Cmpt, PatchField, GeoMesh>& yy,
    GeometricField<Cmpt, PatchField, GeoMesh>& yz,
    GeometricField<Cmpt, PatchField, GeoMesh>& zz
)
{
    Foam::unzip
    (
        input.primitiveField(),
        xx.primitiveFieldRef(), xy.primitiveFieldRef(), xz.primitiveFieldRef(),
        yy.primitiveFieldRef(), yz.primitiveFieldRef(),
        zz.primitiveFieldRef()
    );

    Foam::unzip
    (
        input.boundaryField(),
        xx.boundaryFieldRef(), xy.boundaryFieldRef(), xz.boundaryFieldRef(),
        yy.boundaryFieldRef(), yz.boundaryFieldRef(),
        zz.boundaryFieldRef()
    );
}


template<class Cmpt, template<class> class PatchField, class GeoMesh>
void Foam::unzipDiag
(
    const GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& result
)
{
    Foam::unzipDiag(input.primitiveField(), result.primitiveFieldRef());

    Foam::unzipDiag(input.boundaryField(), result.boundaryFieldRef());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(symmTensor, vector, sqr, sqr)
UNARY_FUNCTION(symmTensor, symmTensor, innerSqr, sqr)

UNARY_FUNCTION(scalar, symmTensor, tr, transform)
UNARY_FUNCTION(sphericalTensor, symmTensor, sph, transform)
UNARY_FUNCTION(symmTensor, symmTensor, symm, transform)
UNARY_FUNCTION(symmTensor, symmTensor, twoSymm, transform)
UNARY_FUNCTION(symmTensor, symmTensor, dev, transform)
UNARY_FUNCTION(symmTensor, symmTensor, dev2, transform)
UNARY_FUNCTION(scalar, symmTensor, det, pow3)
UNARY_FUNCTION(symmTensor, symmTensor, cof, pow2)
UNARY_FUNCTION(symmTensor, symmTensor, inv, inv)


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, symmTensor, *, hdual, transform)

BINARY_OPERATOR(tensor, symmTensor, symmTensor, &, '&', dot)
BINARY_TYPE_OPERATOR(tensor, symmTensor, symmTensor, &, '&', dot)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
