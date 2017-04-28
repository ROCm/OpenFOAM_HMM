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

Class
    Foam::Random

Description
    Random number generator.


SourceFiles
    RandomI.H
    Random.C
    RandomTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Random_H
#define Random_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Random;

/*---------------------------------------------------------------------------*\
                           Class Random Declaration
\*---------------------------------------------------------------------------*/

class Random
{
    // Private data

        //- Buffer used by re-entrant random number generator
        List<char> buffer_;

        //- Initial random number seed
        label seed_;

        //- Indicator, which tells if there is a stored gaussian sample
        bool hasGaussSample_;

        //- Stored sample value
        scalar gaussSample_;


    // Private Member Functions

        //- Returns the current sample
        scalar scalar01();


public:

    // Constructors

        //- Construct given seed and sample count
        Random(const label seed = 123456);

        //- Copy constructor with optional reset of sampleI
        Random(const Random& r, const bool reset = false);


    // Destructor
    ~Random();


    // Member functions

        // Access

            //- Return const access to the initial random number seed
            inline label seed() const;

            //- Reset the random number generator seed
            void reset(const label seed);


        // Evaluation

            // Random numbers

                //- Return a random bit
                int bit();

                //- Return a sample whose components lie in the range 0-1
                template<class Type>
                Type sample01();

                //- Return a sample whose components are normally distributed
                //  with zero mean and unity variance N(0, 1)
                template<class Type>
                Type GaussNormal();

                //- Return a sample between start and end
                template<class Type>
                Type position(const Type& start, const Type& end);

                //- Randomise value in the range 0-1
                template<class Type>
                void randomise01(Type& value);

                //- Shuffle the values in the list
                template<class Type>
                void shuffle(UList<Type>& values);


            // Global random numbers - consistent across all processors

                //- Return a sample whose components lie in the range 0-1
                template<class Type>
                Type globalSample01();

                //- Return a sample whose components are normally distributed
                //  with zero mean and unity variance N(0, 1)
                template<class Type>
                Type globalGaussNormal();

                //- Return a sample between start and end
                template<class Type>
                Type globalPosition(const Type& start, const Type& end);

                //- Randomise value in the range 0-1
                template<class Type>
                void globalRandomise01(Type& value);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Template specialisations

template<>
scalar Random::sample01<scalar>();

template<>
label Random::sample01<label>();

template<>
scalar Random::GaussNormal<scalar>();

template<>
label Random::GaussNormal<label>();

template<>
scalar Random::position<scalar>
(
    const scalar& start,
    const scalar& end
);

template<>
label Random::position<label>(const label& start, const label& end);

template<>
scalar Random::globalSample01<scalar>();

template<>
label Random::globalSample01<label>();

template<>
scalar Random::globalGaussNormal<scalar>();

template<>
label Random::globalGaussNormal<label>();

template<>
scalar Random::globalPosition<scalar>
(
    const scalar& start,
    const scalar& end
);

template<>
label Random::globalPosition<label>(const label& start, const label& end);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RandomI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "RandomTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
