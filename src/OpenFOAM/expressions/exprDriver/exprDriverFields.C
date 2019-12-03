/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "exprDriver.H"
#include "Tuple2.H"
#include "FieldOps.H"
#include "Random.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//! \cond file-scope

template<class T1, class T2>
static Foam::Tuple2<T1,T2> findMinData
(
    const Field<T1>& vals,
    const Field<T2>& data
)
{
    typedef Tuple2<T1,T2> retType;

    retType result(pTraits<T1>::max, Zero);

    const label i = findMin(vals);
    if (i != -1)
    {
        result.first() = vals[i];
        result.second() = data[i];
    }

    Foam::combineReduce(result, minFirstEqOp<T1>());
    return result;
}


template<class T1, class T2>
static Foam::Tuple2<T1,T2> findMaxData
(
     const Field<T1>& vals,
     const Field<T2>& data
)
{
    typedef Tuple2<T1,T2> retType;

    retType result(pTraits<T1>::min, Zero);

    const label i = findMax(vals);
    if (i != -1)
    {
        result.first() = vals[i];
        result.second() = data[i];
    }

    Foam::combineReduce(result, maxFirstEqOp<T1>());
    return result;
}

//! \endcond
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::expressions::exprDriver::fill_random
(
    scalarField& field,
    label seed,
    const bool gaussian
) const
{
    if (gaussian)
    {
        Random::gaussianGeneratorOp<scalar> gen(seed);

        FieldOps::assign(field, field, gen);
    }
    else
    {
        Random::uniformGeneratorOp<scalar> gen(seed);

        FieldOps::assign(field, field, gen);
    }
}


Foam::point Foam::expressions::exprDriver::getPositionOfMinimum
(
    const scalarField& vals,
    const pointField& locs
)
{
    return findMinData(vals, locs).second();
}


Foam::point Foam::expressions::exprDriver::getPositionOfMaximum
(
    const scalarField& vals,
    const pointField& locs
)
{
    return findMaxData(vals, locs).second();
}


// ************************************************************************* //
