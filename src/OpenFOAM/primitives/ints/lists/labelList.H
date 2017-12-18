/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::labelList

Description
    A List of labels.

Typedef
    Foam::labelSubList

Description
    A SubList of labels.

Typedef
    Foam::labelUList

Description
    A UList of labels.

Typedef
    Foam::labelListList

Description
    A List of labelList.

Typedef
    Foam::labelListListList

Description
    A List of labelListList.

\*---------------------------------------------------------------------------*/

#ifndef labelList_H
#define labelList_H

#include "label.H"
#include "List.H"
#include "SubList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // NB: labelUList is defined in UList itself
    typedef List<label> labelList;
    typedef SubList<label> labelSubList;
    typedef List<labelList> labelListList;
    typedef List<labelListList> labelListListList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
