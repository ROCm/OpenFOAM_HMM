/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, licensed under GNU General Public License
    <http://www.gnu.org/licenses/>.

Typedef
    Foam::meshSubsetHelper

Description
    Compatibility name. Rename (JUL-2018) as Foam::fvMeshSubsetProxy

\*---------------------------------------------------------------------------*/

#ifndef meshSubsetHelper_H
#define meshSubsetHelper_H

#include "fvMeshSubsetProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
     typedef fvMeshSubsetProxy meshSubsetHelper;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
