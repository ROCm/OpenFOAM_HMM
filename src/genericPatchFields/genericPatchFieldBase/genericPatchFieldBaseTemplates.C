/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "genericPatchFieldBase.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class MapperType>
void Foam::genericPatchFieldBase::mapGeneric
(
    const genericPatchFieldBase& rhs,
    const MapperType& mapper
)
{
    forAllConstIters(rhs.scalarFields_, iter)
    {
        scalarFields_.insert
        (
            iter.key(),
            autoPtr<scalarField>::New(*iter(), mapper)
        );
    }

    forAllConstIters(rhs.vectorFields_, iter)
    {
        vectorFields_.insert
        (
            iter.key(),
            autoPtr<vectorField>::New(*iter(), mapper)
        );
    }

    forAllConstIters(rhs.sphTensorFields_, iter)
    {
        sphTensorFields_.insert
        (
            iter.key(),
            autoPtr<sphericalTensorField>::New(*iter(), mapper)
        );
    }

    forAllConstIters(rhs.symmTensorFields_, iter)
    {
        symmTensorFields_.insert
        (
            iter.key(),
            autoPtr<symmTensorField>::New(*iter(), mapper)
        );
    }

    forAllConstIters(rhs.tensorFields_, iter)
    {
        tensorFields_.insert
        (
            iter.key(),
            autoPtr<tensorField>::New(*iter(), mapper)
        );
    }
}


template<class MapperType>
void Foam::genericPatchFieldBase::autoMapGeneric
(
    const MapperType& mapper
)
{
    forAllIters(scalarFields_, iter)
    {
        (*iter)->autoMap(mapper);
    }

    forAllIters(vectorFields_, iter)
    {
        (*iter)->autoMap(mapper);
    }

    forAllIters(sphTensorFields_, iter)
    {
        (*iter)->autoMap(mapper);
    }

    forAllIters(symmTensorFields_, iter)
    {
        (*iter)->autoMap(mapper);
    }

    forAllIters(tensorFields_, iter)
    {
        (*iter)->autoMap(mapper);
    }
}


// ************************************************************************* //
