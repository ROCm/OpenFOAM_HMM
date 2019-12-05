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
    This file is part of OpenFOAM, licensed under GNU General Public License
    <http://www.gnu.org/licenses/>.

Typedef
    Foam::PackedBoolList

Description
    Compatibility name. Superseded (MAR-2018) by Foam::bitSet

\*---------------------------------------------------------------------------*/

#ifndef PackedBoolList_H
#define PackedBoolList_H

#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
     typedef bitSet PackedBoolList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
