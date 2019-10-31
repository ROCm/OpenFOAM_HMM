/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::labelPairPair

Description
    A pair of labelPairs.

Typedef
    Foam::labelPairList

Description
    List of labelPairs.

Typedef
    Foam::labelPairUList

Description
    UList of labelPairs.

Typedef
    Foam::labelPairSubList

Description
    A SubList of labelPairs.

\*---------------------------------------------------------------------------*/

#ifndef labelPair_H
#define labelPair_H

#include "List.H"
#include "SubList.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef Pair<labelPair> labelPairPair;
    typedef List<labelPair> labelPairList;
    typedef UList<labelPair> labelPairUList;
    typedef SubList<labelPair> labelPairSubList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
