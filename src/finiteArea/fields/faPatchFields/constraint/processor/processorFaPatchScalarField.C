/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "processorFaPatchScalarField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::processorFaPatchField<Foam::scalar>::transformCoupleField
(
    solveScalarField& f,
    const direction cmpt
) const
{}


template<>
void Foam::processorFaPatchField<Foam::scalar>::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    procPatch_.send
    (
        commsType,
        patch().patchInternalField(psiInternal)()
    );
}


template<>
void Foam::processorFaPatchField<Foam::scalar>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField&,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    solveScalarField pnf
    (
        procPatch_.receive<solveScalar>(commsType, this->size())()
    );

    const labelUList& edgeFaces = patch().edgeFaces();

    if (add)
    {
        forAll(edgeFaces, facei)
        {
            result[edgeFaces[facei]] -= coeffs[facei]*pnf[facei];
        }
    }
    else
    {
        forAll(edgeFaces, facei)
        {
            result[edgeFaces[facei]] -= coeffs[facei]*pnf[facei];
        }
    }
}


// ************************************************************************* //
