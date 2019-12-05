/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::Instant

Description
    A tuple of value and key.
    The value often corresponds to a time value, thus the naming of the class.
    The key will usually be a time name or a file name etc.

SourceFiles
    Instant.C

\*---------------------------------------------------------------------------*/

#ifndef Instant_H
#define Instant_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Instant Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Instant
{
    // Private Data

        //- The value (eg, time)
        scalar val_;

        //- The name/key
        T key_;


public:

    // Public classes

        //- Less function for sorting
        struct less
        {
            bool operator()(const Instant& a, const Instant& b) const
            {
                return a.value() < b.value();
            }
        };


    // Constructors

        //- Construct null, with value = 0.
        Instant();

        //- Copy construct from components
        Instant(scalar val, const T& key);

        //- Move construct from components
        Instant(scalar val, T&& key);

        //- Copy construct
        Instant(const Instant&) = default;

        //- Move construct
        Instant(Instant&&) = default;


    // Member Functions

        //- The value (const access)
        scalar value() const
        {
            return val_;
        }

        //- The value (non-const access)
        scalar& value()
        {
            return val_;
        }

        //- The name/key (const access)
        const T& name() const
        {
            return key_;
        }

        //- The name/key (non-const access)
        T& name()
        {
            return key_;
        }

        //- True if values are equal (includes SMALL for rounding)
        bool equal(scalar val) const;

        //- True if values are equal (includes SMALL for rounding)
        template<class T2>
        bool equal(const Instant<T2>& other) const;


    // Operators

        //- Copy assignment
        Instant& operator=(const Instant&) = default;

        //- Move assignment
        Instant& operator=(Instant&&) = default;
};


// IOstream Operators

template<class T> Istream& operator>>(Istream& is, Instant<T>& inst);
template<class T> Ostream& operator<<(Ostream& os, const Instant<T>& inst);


// Global Operators

//- Compare instant values for equality
template<class T1, class T2>
bool operator==(const Instant<T1>& a, const Instant<T2>& b);

//- Compare instant values for inequality
template<class T1, class T2>
bool operator!=(const Instant<T1>& a, const Instant<T2>& b);

//- Compare instant values for less-than
template<class T1, class T2>
bool operator<(const Instant<T1>& a, const Instant<T2>& b);

//- Compare instant values for greater-than
template<class T1, class T2>
bool operator>(const Instant<T1>& a, const Instant<T2>& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Instant.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
