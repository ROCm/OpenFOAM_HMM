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

#include "LduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
bool Foam::LduMatrix<Type, DType, LUType>::solverPerformance::singular
(
    const Type& wApA
)
{
    for(direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        singular_[cmpt] = component(wApA, cmpt) < vsmall_;
    }

    return singular();
}


template<class Type, class DType, class LUType>
bool Foam::LduMatrix<Type, DType, LUType>::solverPerformance::singular() const
{
    for(direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        if (!singular_[cmpt]) return false;
    }

    return true;
}


template<class Type, class DType, class LUType>
bool Foam::LduMatrix<Type, DType, LUType>::solverPerformance::converged
(
    const Type& Tolerance,
    const Type& RelTolerance
)
{
    if (debug >= 2)
    {
        Info<< solverName_
            << ":  Iteration " << noIterations_
            << " residual = " << finalResidual_
            << endl;
    }

    if
    (
        finalResidual_ < Tolerance
     || (
            RelTolerance > small_*pTraits<Type>::one
         && finalResidual_ < cmptMultiply(RelTolerance, initialResidual_)
        )
    )
    {
        converged_ = true;
    }
    else
    {
        converged_ = false;
    }

    return converged_;
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::solverPerformance::print
(
    Ostream& os
) const
{
    for(direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        os  << solverName_ << ":  Solving for "
            << word(fieldName_ + pTraits<Type>::componentNames[cmpt]);

        if (singular_[cmpt])
        {
            os  << ":  solution singularity" << endl;
        }
        else
        {
            os  << ", Initial residual = " << component(initialResidual_, cmpt)
                << ", Final residual = " << component(finalResidual_, cmpt)
                << ", No Iterations " << noIterations_
                << endl;
        }
    }
}


// ************************************************************************* //
