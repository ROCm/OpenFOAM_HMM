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
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Typedef
    Foam::tmpNrc

Description
    Compatibility name. Superseded (JUL-2020) by Foam::refPtr

\*---------------------------------------------------------------------------*/

#ifndef tmpNrc_H
#define tmpNrc_H

#include "refPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef refPtr tmpNrc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
