/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::functionEntries::evalEntry

Description
    Uses stringOps::toScalar to evaluate mathematical expressions.

    The input can any form of string or, for convenience,
    a '{}' delimited string literal.  In all cases, C/C++ comment stripping
    is also performed.

    For example,

    \verbatim
    a 1;
    b 3;
    c #eval "sin(pi()*$a/$b)";

    d #eval{
        // ignore: sin(pi()*$a/$b)
        sin(degToRad(45))
    };
    \endverbatim

SourceFiles
    evalEntry.C

\*---------------------------------------------------------------------------*/

#ifndef evalEntry_H
#define evalEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                          Class evalEntry Declaration
\*---------------------------------------------------------------------------*/

class evalEntry
:
    public functionEntry
{

    //- Evaluate and return a scalar
    static scalar evaluate(const dictionary& parentDict, Istream& is);


public:

    //- Execute in a primitiveEntry context
    static bool execute
    (
        const dictionary& parentDict,
        primitiveEntry& thisEntry,
        Istream& is
    );

    //- Execute in a sub-dict context
    static bool execute(dictionary& parentDict, Istream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
