/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvPatch::patchInternalField
(
    const UList<Type>& f,
    const labelUList& faceCells,
    Field<Type>& pfld
) const
{
    pfld.resize(size());

    forAll(pfld, i)
    {
        pfld[i] = f[faceCells[i]];
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvPatch::patchInternalField
(
    const UList<Type>& f
) const
{
    auto tpfld = tmp<Field<Type>>::New(size());
    patchInternalField(f, this->faceCells(), tpfld.ref());
    return tpfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvPatch::patchInternalField
(
    const UList<Type>& f,
    const labelUList& faceCells
) const
{
    auto tpfld = tmp<Field<Type>>::New(size());
    patchInternalField(f, faceCells, tpfld.ref());
    return tpfld;
}


template<class Type>
void Foam::fvPatch::patchInternalField
(
    const UList<Type>& f,
    Field<Type>& pfld
) const
{
    patchInternalField(f, this->faceCells(), pfld);
}


template<class GeometricField, class AnyType>
const typename GeometricField::Patch& Foam::fvPatch::patchField
(
    const GeometricField& gf
) const
{
    return gf.boundaryField()[this->index()];
}


// ************************************************************************* //
