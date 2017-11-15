/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::Pair

Description
    An ordered pair of two objects of type \<T\> with first() and second()
    elements.

SourceFiles
    PairI.H

See also
    Foam::Tuple2 for storing two objects of dissimilar types.

\*---------------------------------------------------------------------------*/

#ifndef Pair_H
#define Pair_H

#include "FixedList.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Pair Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Pair
:
    public FixedList<T, 2>
{

public:

    // Constructors

        //- Null constructor
        inline Pair();

        //- Construct from components
        inline Pair(const T& f, const T& s);

        //- Construct from FixedList
        inline Pair(const FixedList<T, 2>& lst);

        //- Construct, optionally sorted with first less-than second
        inline Pair(const T& f, const T& s, const bool doSort);

        //- Construct, optionally sorted with first less-than second
        inline Pair(const FixedList<T, 2>& lst, const bool doSort);

        //- Construct from Istream
        inline Pair(Istream& is);


    // Member Functions

      // Access

        //- Return first element
        using FixedList<T, 2>::first;

        //- Return last element
        using FixedList<T, 2>::last;

        //- Return second element, which is also the last element
        inline const T& second() const;

        //- Return second element, which is also the last element
        inline T& second();

        //- Return other element
        inline const T& other(const T& a) const;


     // Queries

        //- True if first() is less-than second()
        inline bool sorted() const;


     // Editing

        //- Flip the Pair in-place.
        inline void flip();

        //- Sort so that first() is less-than second()
        inline void sort();


     // Comparison

        //- Compare Pairs
        //  \return
        //  -  0: different
        //  - +1: identical values and order used
        //  - -1: identical values, but in reversed order
        static inline int compare(const Pair<T>& a, const Pair<T>& b);
};


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

//- Return reverse of a Pair
template<class T>
Pair<T> reverse(const Pair<T>& p)
{
    return Pair<T>(p.second(), p.first());
}


template<class T>
bool operator==(const Pair<T>& a, const Pair<T>& b)
{
    return (a.first() == b.first() && a.second() == b.second());
}


template<class T>
bool operator!=(const Pair<T>& a, const Pair<T>& b)
{
    return !(a == b);
}


template<class T>
bool operator<(const Pair<T>& a, const Pair<T>& b)
{
    return
    (
        a.first() < b.first()
     ||
        (
            !(b.first() < a.first())
         && a.second() < b.second()
        )
    );
}


template<class T>
bool operator<=(const Pair<T>& a, const Pair<T>& b)
{
    return !(b < a);
}


template<class T>
bool operator>(const Pair<T>& a, const Pair<T>& b)
{
    return (b < a);
}


template<class T>
bool operator>=(const Pair<T>& a, const Pair<T>& b)
{
    return !(a < b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PairI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
