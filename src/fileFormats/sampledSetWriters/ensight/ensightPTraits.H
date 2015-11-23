/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

Class
    Foam::ensightPTraits

Description
    Conversion of OpenFOAM pTraits into the Ensight equivalent

\*---------------------------------------------------------------------------*/

#ifndef ensightPTraits_H
#define ensightPTraits_H

#include "pTraits.H"
#include "fieldTypes.H"
#include "direction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ensightPTraits Declaration
\*---------------------------------------------------------------------------*/

template<class PrimitiveType>
class ensightPTraits
{
public:

    // Static data members

        static const char* const typeName;

        //- Ordering table: return OpenFOAM component given Ensight component
        //  This is used for the symmTensor ordering: Ensight requires
        //  xx yy zz xy xz yz
        static const direction componentOrder[];
};

template<>
const char* const ensightPTraits<scalar>::typeName;

template<>
const direction ensightPTraits<scalar>::componentOrder[];

template<>
const char* const ensightPTraits<vector>::typeName;

template<>
const direction ensightPTraits<vector>::componentOrder[];

template<>
const char* const ensightPTraits<sphericalTensor>::typeName;

template<>
const direction ensightPTraits<sphericalTensor>::componentOrder[];

template<>
const char* const ensightPTraits<symmTensor>::typeName;

template<>
const direction ensightPTraits<symmTensor>::componentOrder[];

template<>
const char* const ensightPTraits<tensor>::typeName;

template<>
const direction ensightPTraits<tensor>::componentOrder[];


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
