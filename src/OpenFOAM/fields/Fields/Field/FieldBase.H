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
    Foam::FieldBase

Description
    Template invariant parts for Field and SubField

SourceFiles
    FieldBase.C

\*---------------------------------------------------------------------------*/

#ifndef FieldBase_H
#define FieldBase_H

#include "refCount.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class FieldBase Declaration
\*---------------------------------------------------------------------------*/

class FieldBase
:
    public refCount
{
public:

    // Static Data Members

        //- Typename for Field
        static const char* const typeName;

        //- Permit read construct from a larger size.
        //  Mostly required for things like column mesh, for example.
        static bool allowConstructFromLargerSize;


    // Constructors

        //- Default construct, refCount zero
        constexpr FieldBase() noexcept
        :
            refCount()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
