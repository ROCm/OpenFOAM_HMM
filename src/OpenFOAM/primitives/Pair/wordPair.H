/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::wordPairPair

Description
    A pair of wordPairs.

Typedef
    Foam::wordPairList

Description
    List of wordPairs.

Typedef
    Foam::wordPairUList

Description
    UList of wordPairs.

Typedef
    Foam::wordPairSubList

Description
    A SubList of wordPairs.

\*---------------------------------------------------------------------------*/

#ifndef wordPair_H
#define wordPair_H

#include "List.H"
#include "SubList.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef Pair<wordPair> wordPairPair;
    typedef List<wordPair> wordPairList;
    typedef UList<wordPair> wordPairUList;
    typedef SubList<wordPair> wordPairSubList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
