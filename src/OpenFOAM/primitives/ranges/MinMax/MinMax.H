/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::MinMax

Description
    A min/max value pair with additional methods.
    In addition to conveniently storing values, it can be used for logic
    operations or to modify data. A few global functions and functors are
    also provided.

    Examples of use.

    Determine min/max limits from a List of values:
    \verbatim
    List<scalar> values = ...;

    // on construction
    MinMax<scalar> range(values);

    range.clear();

    range += val;

    // global minMax() function
    Info<< minMax(values) << nl;
    \endverbatim

    General comparison operations
    \verbatim
    scalar val;
    if (val < range) ... value is below range min
    if (range.contains(val)) ... value within range
    if (range.compare(val) > 0) ... value is above range max
    if (range(val)) ... value within range - as predicate
    \endverbatim

    Since the range has a predicate form, it can be used as a filter method.
    For example,
    \verbatim
    Info<< "values in range: " << subsetList(values, range) << nl;

    boolList mask = ListOps::create<bool>(values, range);
    Info<< "values values in range " << mask << nl;
    \endverbatim

    One particular advantage offered by MinMax is to clip or limit values
    to a particular range. For example,
    \verbatim
    scalarMinMax range(lower, upper);

    scalar val;
    val = range.clip(val)    .. return clip values

    // vs.
    val = min(max(value, lower, upper))
    \endverbatim

    Or when working on lists, the values can be limited in a single pass
    of the data without intermediate memory allocation.
    \verbatim
    scalarField values = ...;

    for (scalar& val : values)
    {
        range.inplaceClip(val);
    }

    // vs.
    values = min(max(values, lower, upper))
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef MinMax_H
#define MinMax_H

#include "scalar.H"
#include "Pair.H"
#include "Tuple2.H"
#include "VectorSpace.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
template<class T> class MinMax;
class zero;

// Common min/max types
typedef MinMax<label> labelMinMax;      //!< A label min/max range
typedef MinMax<scalar> scalarMinMax;    //!< A scalar min/max range


/*---------------------------------------------------------------------------*\
                         Class MinMax Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class MinMax
:
    public Tuple2<T,T>
{
public:

    // Typedefs

        //- The value type the MinMax represents
        typedef T value_type;

        //- Component type
        typedef typename pTraits<T>::cmptType cmptType;


    // Constructors

        //- Construct inverted range
        inline MinMax();

        //- Copy construct from components
        inline MinMax(const T& minVal, const T& maxVal);

        //- Copy construct from components
        inline MinMax(const std::pair<T,T>& range);

        //- Copy construct from components
        inline MinMax(const Pair<T>& range);

        //- Construct with a single zero value
        inline explicit MinMax(const Foam::zero);

        //- Construct with a single initial value
        inline explicit MinMax(const T& val);

        //- Construct from list of values
        inline explicit MinMax(const UList<T>& vals);


    // Static Member Functions

        //- A semi-infinite range from minVal to the type max
        inline static MinMax<T> ge(const T& minVal);

        //- A semi-infinite range from type min to maxVal
        inline static MinMax<T> le(const T& maxVal);

        //- A 0-1 range corresponding to the pTraits zero, one
        inline static MinMax<T> zero_one();


    // Member Functions

    // Access

        //- The min value (first)
        inline const T& min() const noexcept;

        //- The min value (first)
        inline T& min() noexcept;

        //- The max value (second)
        inline const T& max() const noexcept;

        //- The max value (second)
        inline T& max() noexcept;

        //- The min/max average value
        inline T centre() const;

        //- The min to max span. Zero if the range is invalid.
        inline T span() const;

        //- The magnitude of the min to max span. Zero if the range is invalid.
        inline scalar mag() const;

        //- Range is empty if it is inverted
        inline bool empty() const;

        //- Range is valid if it is not inverted
        inline bool valid() const;

        //- Reset to an invalid, inverted range
        inline void clear();


    // Testing / Query

        //- Intersect (union) with the second range.
        //  \return True if the resulting intersection is non-empty.
        inline bool intersect(const MinMax<T>& b);

        //- Test if the ranges overlap
        inline bool overlaps(const MinMax<T>& b) const;

        //- Compares the min/max range with the specified value.
        //  \return
        //  -  0: value is within the range, or range is invalid
        //  - -1: range (max) is less than the value
        //  - +1: range (min) is greater than value
        inline int compare(const T& val) const;

        //- True if the value is within the range
        inline bool contains(const T& val) const;

        //- If out of range, return the respective min/max limits, otherwise
        //- return the value itself.
        //  If the range is invalid, always return the value.
        inline const T& clip(const T& val) const;

        //- Inplace clip value by the min/max limits
        //  \return True if clipping was applied.
        inline bool inplaceClip(T& val) const;


    // Manipulate

        //- Extend the range to include the other min/max range
        inline MinMax<T>& add(const MinMax& other);

        //- Include the value into the range
        inline MinMax<T>& add(const T& val);

        //- Include the values into the range
        inline MinMax<T>& add(const UList<T>& vals);


    // Member Operators

        //- Identical to contains(), for use as a predicate.
        inline bool operator()(const T& val) const;

        //- Extend min/max range to include other range
        //  Can be used in a reduction operation.
        inline MinMax<T>& operator+=(const MinMax<T>& b);

        //- Extend min/max range to include value
        inline MinMax<T>& operator+=(const T& val);

        //- Extend min/max range to include all values
        inline MinMax<T>& operator+=(const UList<T>& vals);

        //- Multiply range by scalar factor
        inline MinMax<T>& operator*=(const scalar& s);

        //- Divide range by scalar factor
        inline MinMax<T>& operator/=(const scalar& s);
};


// Global Functions

//- Min/max range as a string
template<class T>
word name(const MinMax<T>& range)
{
    return '(' + Foam::name(range.min()) + ',' + Foam::name(range.max()) + ')';
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MinMaxI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// Global Functions and Operators
#include "MinMaxOps.H"


// ************************************************************************* //
