/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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



#ifndef OMP_UNIFIED_MEMORY_REQUIRED
#pragma omp requires unified_shared_memory
#define OMP_UNIFIED_MEMORY_REQUIRED
#endif

#include "AtomicAccumulator.H"
#include "macros.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::lduInterfaceField::addToInternalField
(
    Field<Type>& result,
    const bool add,
    const labelUList& faceCells,
    const scalarField& coeffs,
    const Field<Type>& vals
) const
{
    

    if (add)
    {
        //forAll(faceCells, elemi)
	const label loop_len = faceCells.size();
        #pragma omp target teams distribute parallel for if(loop_len>5000)
        for (label elemi = 0; elemi < loop_len; ++elemi)
        {
            atomicAccumulator(result[faceCells[elemi]]) += (coeffs[elemi]*vals[elemi]);
        }
    }
    else
    {
        //forAll(faceCells, elemi)
        const label loop_len = faceCells.size();
        #pragma omp target teams distribute parallel for if(loop_len>5000)
        for (label elemi = 0; elemi < loop_len; ++elemi)
        {
            atomicAccumulator(result[faceCells[elemi]]) -= (coeffs[elemi]*vals[elemi]);
        }
    }
}


// ************************************************************************* //
