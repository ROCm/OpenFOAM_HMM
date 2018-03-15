/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

Description
    Helpers for OpenFOAM reader interfaces in ParaView.

\*---------------------------------------------------------------------------*/

#ifndef foamPvFields_H
#define foamPvFields_H

#include "symmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class foamPvFields Declaration
\*---------------------------------------------------------------------------*/

struct foamPvFields
{
    //- Remapping for some data types
    template<class Type>
    inline static void remapTuple(float vec[])
    {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Template specialization for symmTensor
template<>
inline void Foam::foamPvFields::remapTuple<Foam::symmTensor>(float vec[])
{
    Foam::Swap(vec[1], vec[3]);   // swap XY <-> YY
    Foam::Swap(vec[2], vec[5]);   // swap XZ <-> ZZ
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
