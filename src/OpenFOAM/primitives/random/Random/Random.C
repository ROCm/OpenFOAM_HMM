/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
#include "OSspecific.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::Random::scalar01()
{
    return osRandomDouble(buffer_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Random::Random(const label seed)
:
    buffer_(osRandomBufferSize()),
    seed_(seed),
    hasGaussSample_(false),
    gaussSample_(0)
{
    // Initialise the random number generator
    osRandomSeed(seed_, buffer_);
}


Foam::Random::Random(const Random& r, const bool reset)
:
    buffer_(r.buffer_),
    seed_(r.seed_),
    hasGaussSample_(r.hasGaussSample_),
    gaussSample_(r.gaussSample_)
{
    if (reset)
    {
        hasGaussSample_ = false;
        gaussSample_ = 0;

        // Re-initialise the samples
        osRandomSeed(seed_, buffer_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Random::~Random()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Random::reset(const label seed)
{
    seed_ = seed;
    osRandomSeed(seed_, buffer_);
}


int Foam::Random::bit()
{
    return osRandomInteger(buffer_) & 1;
}


template<>
Foam::scalar Foam::Random::sample01()
{
    return scalar01();
}


template<>
Foam::label Foam::Random::sample01()
{
    return round(scalar01());
}


template<>
Foam::scalar Foam::Random::GaussNormal()
{
    if (hasGaussSample_)
    {
        hasGaussSample_ = false;
        return gaussSample_;
    }
    else
    {
        scalar rsq, v1, v2;
        do
        {
            v1 = 2*scalar01() - 1;
            v2 = 2*scalar01() - 1;
            rsq = sqr(v1) + sqr(v2);
        } while (rsq >= 1 || rsq == 0);

        scalar fac = sqrt(-2*log(rsq)/rsq);

        gaussSample_ = v1*fac;
        hasGaussSample_ = true;

        return v2*fac;
    }
}


template<>
Foam::label Foam::Random::GaussNormal()
{
    return round(GaussNormal<scalar>());
}


template<>
Foam::scalar Foam::Random::position
(
    const scalar& start,
    const scalar& end
)
{
    return start + scalar01()*(end - start);
}


template<>
Foam::label Foam::Random::position(const label& start, const label& end)
{
    return start + round(scalar01()*(end - start));
}


template<>
Foam::scalar Foam::Random::globalSample01()
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = scalar01();
    }

    Pstream::scatter(value);

    return value;
}


template<>
Foam::label Foam::Random::globalSample01()
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = scalar01();
    }

    Pstream::scatter(value);

    return round(value);
}


template<>
Foam::scalar Foam::Random::globalGaussNormal()
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = GaussNormal<scalar>();
    }

    Pstream::scatter(value);

    return value;
}


template<>
Foam::label Foam::Random::globalGaussNormal()
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = GaussNormal<scalar>();
    }

    Pstream::scatter(value);

    return round(value);
}


template<>
Foam::scalar Foam::Random::globalPosition
(
    const scalar& start,
    const scalar& end
)
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = scalar01()*(end - start);
    }

    Pstream::scatter(value);

    return start + value;
}


template<>
Foam::label Foam::Random::globalPosition
(
    const label& start,
    const label& end
)
{
    label value = labelMin;

    if (Pstream::master())
    {
        value = round(scalar01()*(end - start));
    }

    Pstream::scatter(value);

    return start + value;
}


// ************************************************************************* //
