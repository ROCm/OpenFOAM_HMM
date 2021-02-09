/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Various functors for unary and binary operations.
    Can be used for parallel combine-reduce operations or other places
    requiring a functor.

\*---------------------------------------------------------------------------*/

#ifndef ops_H
#define ops_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Assignment operation taking two parameters, returning void.
// Alters the value of the first parameter.
//     Eg, plusEqOp for (x += y)
#define EqOp(opName, op)                                                       \
                                                                               \
    template<class T1, class T2>                                               \
    struct opName##Op2                                                         \
    {                                                                          \
        void operator()(T1& x, const T2& y) const                              \
        {                                                                      \
            op;                                                                \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T>                                                          \
    struct opName##Op                                                          \
    {                                                                          \
        void operator()(T& x, const T& y) const                                \
        {                                                                      \
            op;                                                                \
        }                                                                      \
    };


EqOp(eq, x = y)
EqOp(plusEq, x += y)
EqOp(minusEq, x -= y)
EqOp(multiplyEq, x *= y)
EqOp(divideEq, x /= y)
EqOp(eqSqr, x = sqr(y))
EqOp(eqMag, x = mag(y))
EqOp(eqMagSqr, x = magSqr(y))
EqOp(plusEqMagSqr, x += magSqr(y))
EqOp(maxEq, x = max(x, y))
EqOp(minEq, x = min(x, y))
EqOp(minMagSqrEq, x = (magSqr(x) <= magSqr(y) ? x : y))
EqOp(maxMagSqrEq, x = (magSqr(x) >= magSqr(y) ? x : y))

EqOp(andEq, x = (x && y))
EqOp(orEq,  x = (x || y))
EqOp(xorEq, x = (x != y))
EqOp(bitAndEq, x = (x & y))
EqOp(bitOrEq,  x = (x | y))
EqOp(bitXorEq, x = (x ^ y))

EqOp(eqMinus, x = -y)

EqOp(nopEq, (void)x)

#undef EqOp


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Warning about unused result
#if __GNUC__
#define WARNRETURN __attribute__((warn_unused_result))
#else
#define WARNRETURN
#endif

// Operation taking two parameters, returning the first type.
// Neither parameter is altered.
//     Eg, plusOp for (x + y)

#define Op(opName, op)                                                         \
                                                                               \
    template<class T, class T1, class T2>                                      \
    struct opName##Op3                                                         \
    {                                                                          \
        T operator()(const T1& x, const T2& y) const WARNRETURN                \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T1, class T2>                                               \
    struct opName##Op2                                                         \
    {                                                                          \
        T1 operator()(const T1& x, const T2& y) const WARNRETURN               \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T>                                                          \
    struct opName##Op                                                          \
    {                                                                          \
        T operator()(const T& x, const T& y) const WARNRETURN                  \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };


// Operations taking two parameters (unaltered), returning bool

#define BoolOp(opName, op)                                                     \
                                                                               \
    template<class T1, class T2>                                               \
    struct opName##Op2                                                         \
    {                                                                          \
        bool operator()(const T1& x, const T2& y) const WARNRETURN             \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };                                                                         \
                                                                               \
    template<class T>                                                          \
    struct opName##Op                                                          \
    {                                                                          \
        bool operator()(const T& x, const T& y) const WARNRETURN               \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };


// Operations taking one parameter, returning bool.
// The comparison value is defined during construction

#define Bool1Op(opName, op)                                                    \
                                                                               \
    template<class T>                                                          \
    struct opName##Op1                                                         \
    {                                                                          \
        const T& value;                                                        \
                                                                               \
        opName##Op1(const T& v) : value(v) {}                                  \
                                                                               \
        bool operator()(const T& x) const WARNRETURN                           \
        {                                                                      \
            return op;                                                         \
        }                                                                      \
    };


// Weighting operations

#define WeightedOp(opName, op)                                                 \
                                                                               \
    template<class T, class CombineOp>                                         \
    class opName##WeightedOp                                                   \
    {                                                                          \
        const CombineOp& cop_;                                                 \
                                                                               \
    public:                                                                    \
                                                                               \
        opName##WeightedOp(const CombineOp& cop)                               \
        :                                                                      \
            cop_(cop)                                                          \
        {}                                                                     \
                                                                               \
        void operator()                                                        \
        (                                                                      \
            T& x,                                                              \
            const label index,                                                 \
            const T& y,                                                        \
            const scalar weight                                                \
        ) const                                                                \
        {                                                                      \
            cop_(x, op);                                                       \
        }                                                                      \
    };                                                                         \


Op(sum, x + y)

Op(plus, x + y)
Op(minus, x - y)
Op(multiply, x * y)
Op(divide, x / y)
Op(cmptMultiply, cmptMultiply(x, y))
Op(cmptPow, cmptPow(x, y))
Op(cmptDivide, cmptDivide(x, y))
Op(stabilise, stabilise(x, y))
Op(max, max(x, y))
Op(min, min(x, y))
Op(minMagSqr, (magSqr(x) <= magSqr(y) ? x : y))
Op(maxMagSqr, (magSqr(x) >= magSqr(y) ? x : y))
Op(minMod, minMod(x, y))

Op(bitAnd, (x & y))
Op(bitOr,  (x | y))
Op(bitXor, (x ^ y))

BoolOp(and, x && y)
BoolOp(or,  x || y)
BoolOp(xor, (!x) != (!y))   // With forced bool context
BoolOp(equal, x == y)
BoolOp(notEqual, x != y)
BoolOp(less, x < y)
BoolOp(lessEq, x <= y)
BoolOp(greater, x > y)
BoolOp(greaterEq, x >= y)

Bool1Op(equal, x == value)
Bool1Op(notEqual, x != value)
Bool1Op(less, x < value)
Bool1Op(lessEq, x <= value)
Bool1Op(greater, x > value)
Bool1Op(greaterEq, x >= value)

WeightedOp(multiply, (weight*y))

#undef Op
#undef BoolOp
#undef Bool1Op
#undef WeightedOp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Three-way comparison operation of two parameters,
//  similar to the \c <=> operator in C++20.
//
//  \return a negative value for less, a positive value for greater,
//      and zero for equal value.
template<class T>
struct compareOp
{
    int operator()(const T& a, const T& b) const WARNRETURN
    {
        return (a < b) ? -1 : (b < a) ? 1 : 0;
    }
};


//- Deprecated(2020-11) use nameOp (word.H)
//  \deprecated(2020-11) use nameOp
template<class T> struct getNameOp : nameOp<T> {};

//- Deprecated(2020-11) use typeOp (word.H)
//  \deprecated(2020-11) use typeOp
template<class T> struct getTypeOp : typeOp<T> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef WARNRETURN

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
