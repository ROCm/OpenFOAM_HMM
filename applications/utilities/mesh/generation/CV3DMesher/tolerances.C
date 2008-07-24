/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "CV3D.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CV3D::tolerances::tolerances
(
    const dictionary& controlDict,
    const scalar minCellSize,
    const boundBox& bb
)
:
    span
    (
        max(mag(bb.max().x()), mag(bb.min().x()))
      + max(mag(bb.max().y()), mag(bb.min().y()))
      + max(mag(bb.max().z()), mag(bb.min().z()))
    ),
    span2(Foam::sqr(span)),

    minEdgeLen(readScalar(controlDict.lookup("minEdgeLenCoeff"))*minCellSize),
    minEdgeLen2(Foam::sqr(minEdgeLen)),

    ppDist(readScalar(controlDict.lookup("ppDistCoeff"))*minCellSize),
    ppDist2(Foam::sqr(ppDist))
{}


// ************************************************************************* //
