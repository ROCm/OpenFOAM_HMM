/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

    As a special exception, you have permission to link this program with the
    CGAL library and distribute executables, as long as you follow the
    requirements of the GNU GPL in regard to all of the software in the
    executable aside from CGAL.

Class
    Foam::vectorTools

Description
    Functions for analysing the relationships between vectors

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef vectorTools_H
#define vectorTools_H

#include "vector.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace vectorTools Declaration
\*---------------------------------------------------------------------------*/

//- Collection of functions for testing relationships between two vectors.
namespace vectorTools
{

    //- Test if a and b are parallel: a.b = 1
    template <typename T>
    bool areParallel
    (
        const Vector<T>& a,
        const Vector<T>& b,
        const T& tolerance = SMALL
    )
    {
        return (mag(a ^ b) < tolerance) ? true : false;
//        return ( mag( mag(a & b)/(mag(a)*mag(b)) - 1.0 ) < tolerance )
//               ? true
//               : false;
    }

    //- Test if a and b are orthogonal: a.b = 0
    template <typename T>
    bool areOrthogonal
    (
        const Vector<T>& a,
        const Vector<T>& b,
        const T& tolerance = SMALL
    )
    {
        return (mag(a & b) < tolerance) ? true : false;
    }

    //- Test if angle between a and b is acute: a.b > 0
    template <typename T>
    bool areAcute
    (
        const Vector<T>& a,
        const Vector<T>& b
    )
    {
        return ((a & b) > 0) ? true : false;
    }

    //- Test if angle between a and b is obtuse: a.b < 0
    template <typename T>
    bool areObtuse
    (
        const Vector<T>& a,
        const Vector<T>& b
    )
    {
        return ((a & b) < 0) ? true : false;
    }

    //- Calculate angle between a and b in radians
    template <typename T>
    T radAngleBetween
    (
        const Vector<T>& a,
        const Vector<T>& b,
        const T& tolerance = SMALL
    )
    {
        return Foam::acos( (a & b)/(mag(a)*mag(b) + tolerance) );
    }

    //- Calculate angle between a and b in degrees
    template <typename T>
    T degAngleBetween
    (
        const Vector<T>& a,
        const Vector<T>& b,
        const T& tolerance = SMALL
    )
    {
        return Foam::radToDeg(radAngleBetween(a, b, tolerance));
    }

} // End namespace vectorTools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
