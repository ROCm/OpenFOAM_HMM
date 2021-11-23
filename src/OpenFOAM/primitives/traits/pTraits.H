/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::pTraits

Description
    A traits class, which is primarily used for primitives.

    All primitives need a specialised version of this class. The
    specialised version will normally also require a conversion
    method.

\*---------------------------------------------------------------------------*/

#ifndef pTraits_H
#define pTraits_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;

/*---------------------------------------------------------------------------*\
                           Class pTraits Declaration
\*---------------------------------------------------------------------------*/

// The base implementation is a pass-through to the base class.
// Accordingly it inherits all static methods (eg, typeName etc).
template<class Base>
class pTraits
:
    public Base
{
public:

    // Constructors

        //- Copy construct from base class
        explicit pTraits(const Base& obj)
        :
            Base(obj)
        {}

        //- Construct from Istream
        explicit pTraits(Istream& is)
        :
            Base(is)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
