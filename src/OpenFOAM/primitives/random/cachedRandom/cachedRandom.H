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
    Foam::cachedRandom

Description
    Random number generator.


SourceFiles
    cachedRandomI.H
    cachedRandom.C
    cachedRandomTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cachedRandom_H
#define cachedRandom_H

#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class cachedRandom;

/*---------------------------------------------------------------------------*\
                       Class cachedRandom Declaration
\*---------------------------------------------------------------------------*/

class cachedRandom
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
        cachedRandom(const label seed = 123456);

        //- Copy constructor with optional reset of sampleI
        cachedRandom(const cachedRandom& cr, const bool reset = false);


    // Destructor
    ~cachedRandom();


    // Member functions

        // Access

            //- Return const access to the initial random number seed
            inline label seed() const;

            void reset(const label seed);

        // Evaluation

            // Random numbers

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
scalar cachedRandom::sample01<scalar>();

template<>
label cachedRandom::sample01<label>();

template<>
scalar cachedRandom::GaussNormal<scalar>();

template<>
label cachedRandom::GaussNormal<label>();

template<>
scalar cachedRandom::position<scalar>
(
    const scalar& start,
    const scalar& end
);

template<>
label cachedRandom::position<label>(const label& start, const label& end);

template<>
scalar cachedRandom::globalSample01<scalar>();

template<>
label cachedRandom::globalSample01<label>();

template<>
scalar cachedRandom::globalGaussNormal<scalar>();

template<>
label cachedRandom::globalGaussNormal<label>();

template<>
scalar cachedRandom::globalPosition<scalar>
(
    const scalar& start,
    const scalar& end
);

template<>
label cachedRandom::globalPosition<label>(const label& start, const label& end);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cachedRandomI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cachedRandomTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
