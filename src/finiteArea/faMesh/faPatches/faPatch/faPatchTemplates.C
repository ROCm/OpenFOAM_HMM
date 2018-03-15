/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

\*----------------------------------------------------------------------------*/

#include "faPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::faPatch::patchInternalField
(
    const UList<Type>& f
) const
{
    tmp<Field<Type>> tpif (new Field<Type>(size()));
    Field<Type>& pif = tpif.ref();

    const labelUList& edgeFaces = this->edgeFaces();

    forAll(pif, facei)
    {
        pif[facei] = f[edgeFaces[facei]];
    }

    return tpif;
}


template<class GeometricField, class Type>
const typename GeometricField::Patch& Foam::faPatch::patchField
(
    const GeometricField& gf
) const
{
    return gf.boundaryField()[index()];
}


// ************************************************************************* //
