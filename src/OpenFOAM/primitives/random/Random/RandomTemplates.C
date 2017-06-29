/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "Random.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Random::sample01()
{
    Type value;
    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        value.component(cmpt) = scalar01();
    }

    return value;
}


template<class Type>
Type Foam::Random::GaussNormal()
{
    Type value;
    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        value.component(cmpt) = GaussNormal<scalar>();
    }

    return value;
}


template<class Type>
Type Foam::Random::position(const Type& start, const Type& end)
{
    Type value(start);
    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        value.component(cmpt) +=
            scalar01()*(end.component(cmpt) - start.component(cmpt));
    }

    return value;
}


template<class Type>
void Foam::Random::randomise01(Type& value)
{
    value = sample01<Type>();
}


template<class Type>
void Foam::Random::shuffle(UList<Type>& values)
{
    const label nSample = values.size();
    label posI = nSample - 1;

    for (label i = 1; i < nSample; i++)
    {
        label j = position<label>(0, posI);
        Type t = values[j];
        values[j] = values[posI];
        values[posI] = t;
        posI--;
    }
}


template<class Type>
Type Foam::Random::globalSample01()
{
    Type value = -GREAT*pTraits<Type>::one;

    if (Pstream::master())
    {
        value = sample01<Type>();
    }

    Pstream::scatter(value);

    return value;
}


template<class Type>
Type Foam::Random::globalGaussNormal()
{
    Type value = -GREAT*pTraits<Type>::one;

    if (Pstream::master())
    {
        value = GaussNormal<Type>();
    }

    Pstream::scatter(value);

    return value;
}


template<class Type>
Type Foam::Random::globalPosition(const Type& start, const Type& end)
{
    Type value = -GREAT*pTraits<Type>::one;

    if (Pstream::master())
    {
        value = position<Type>(start, end);
    }

    Pstream::scatter(value);

    return value;
}


template<class Type>
void Foam::Random::globalRandomise01(Type& value)
{
    value = -GREAT*pTraits<Type>::one;

    if (Pstream::master())
    {
        value = sample01<Type>();
    }

    Pstream::scatter(value);
}


// ************************************************************************* //
