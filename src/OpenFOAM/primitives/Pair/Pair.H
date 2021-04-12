/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class Pair;

// Common pair types
typedef Pair<label> labelPair;      //!< A pair of labels
typedef Pair<word> wordPair;        //!< A pair of words


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

        //- Default construct
        Pair() = default;

        //- Copy construct from components
        inline Pair(const T& f, const T& s);

        //- Move construct from components
        inline Pair(T&& f, T&& s);

        //- Copy construct from std::pair
        inline Pair(const std::pair<T,T>& vals);

        //- Move construct from std::pair
        inline Pair(std::pair<T,T>&& vals);

        //- Copy construct FixedList of two items
        inline Pair(const FixedList<T, 2>& list);

        //- Copy construct, optionally sorted with first less-than second
        inline Pair(const T& f, const T& s, const bool doSort);

        //- Copy construct, optionally sorted with first less-than second
        inline Pair(const FixedList<T, 2>& list, const bool doSort);

        //- Construct from Istream
        inline explicit Pair(Istream& is);


    // Member Functions

    // Access

        //- Return first element
        using FixedList<T, 2>::first;

        //- Return last element
        using FixedList<T, 2>::last;

        //- Return second element, which is also the last element
        inline const T& second() const noexcept;

        //- Return second element, which is also the last element
        inline T& second() noexcept;

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


    // Hashing

        //- Symmetric hashing functor for Pair, hashes lower value first
        //  Regular hasher inherited from FixedList
        struct symmHasher
        {
            unsigned operator()(const Pair<T>& obj, unsigned seed=0) const
            {
                Foam::Hash<T> op;
                if (obj.second() < obj.first())
                {
                    return op(obj.first(), op(obj.second(), seed));
                }
                else
                {
                    return op(obj.second(), op(obj.first(), seed));
                }
            }
        };
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Pair is contiguous if the type is contiguous
template<class T>
struct is_contiguous<Pair<T>> : is_contiguous<T> {};

//- Check for Pair of labels
template<class T>
struct is_contiguous_label<Pair<T>> : is_contiguous_label<T> {};

//- Check for Pair of scalars
template<class T>
struct is_contiguous_scalar<Pair<T>> : is_contiguous_scalar<T> {};

//- Hashing for Pair of data
template<class T>
struct Hash<Pair<T>> : Pair<T>::hasher {};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Return reverse of a Pair
template<class T>
Pair<T> reverse(const Pair<T>& p)
{
    return Pair<T>(p.second(), p.first());
}


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

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
     || (!(b.first() < a.first()) && a.second() < b.second())
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
