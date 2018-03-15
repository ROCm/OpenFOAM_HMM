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

Typedef
    Foam::doubleScalar

Description
    Floating-point double precision scalar type.

SourceFiles
    doubleScalar.C

\*---------------------------------------------------------------------------*/

#ifndef doubleScalar_H
#define doubleScalar_H

#include "doubleFloat.H"
#include "direction.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef double doubleScalar;

// Largest and smallest scalar values allowed in certain parts of the code.
// (15 is the number of significant figures in an
//  IEEE double precision number.  See limits.h or float.h)
constexpr doubleScalar doubleScalarGREAT = 1.0e+15;
constexpr doubleScalar doubleScalarVGREAT = 1.0e+300;
constexpr doubleScalar doubleScalarROOTVGREAT = 1.0e+150;
constexpr doubleScalar doubleScalarSMALL = 1.0e-15;
constexpr doubleScalar doubleScalarROOTSMALL = 3.0e-8;
constexpr doubleScalar doubleScalarVSMALL = 1.0e-300;
constexpr doubleScalar doubleScalarROOTVSMALL = 1.0e-150;


#define Scalar doubleScalar
#define ScalarVGREAT doubleScalarVGREAT
#define ScalarVSMALL doubleScalarVSMALL
#define ScalarROOTVGREAT doubleScalarROOTVGREAT
#define ScalarROOTVSMALL doubleScalarROOTVSMALL
#define ScalarRead readDouble


inline Scalar mag(const Scalar s)
{
    return ::fabs(s);
}

inline Scalar hypot(const Scalar x, const Scalar y)
{
    return ::hypot(x, y);
}

inline Scalar atan2(const Scalar y, const Scalar x)
{
    return ::atan2(y, x);
}

inline Scalar jn(const int n, const Scalar s)
{
    return ::jn(n, s);
}

inline Scalar yn(const int n, const Scalar s)
{
    return ::yn(n, s);
}

// Normal (double-precision) transcendental functions
#define transFunc(func)            \
inline Scalar func(const Scalar s) \
{                                  \
    return ::func(s);              \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Scalar.H"

#undef Scalar
#undef ScalarVGREAT
#undef ScalarVSMALL
#undef ScalarROOTVGREAT
#undef ScalarROOTVSMALL
#undef ScalarRead
#undef transFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
