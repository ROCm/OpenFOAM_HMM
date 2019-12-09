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
    Uses expressions::fieldExprDriver to evaluate mathematical expressions
    with scalars, vectors etc.

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

Note
    The string expansions support use of environment variables.
    Unknown variables will expand to an empty string, so it can be advisable
    to an expansion with an alternative. For example,

    \verbatim
    d #eval{ sin(degToRad( ${angle:-0} )) };
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

    //- Evaluate and return a token list
    static tokenList evaluate(const dictionary& parentDict, Istream& is);


public:

    //- Execute in a primitiveEntry context
    static bool execute
    (
        const dictionary& parentDict,
        primitiveEntry& thisEntry,
        Istream& is
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
