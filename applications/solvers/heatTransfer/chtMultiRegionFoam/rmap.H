/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

Global
    rmap

Description
    map field on subset of mesh onto overall field. Rewrites boundary
    conditions to be mixed.

    The source fields can have different patch types for the same destination
    patch. To work around this it attempts to convert all patch fields into
    mixed type since this can accomodate anything from fixedValue to
    fixedGradient.

SourceFiles
    rmap.C

\*---------------------------------------------------------------------------*/

#ifndef rmap_H
#define rmap_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Map patchField
template<class Type>
static void rmap
(
    fvPatchField<Type>& destBC,
    const labelList& reverseAddressing,
    const fvPatchField<Type>& sourceBC
);

//- Map volField
template<class Type>
static void rmap
(
    GeometricField<Type, fvPatchField, volMesh>& dest,
    const GeometricField<Type, fvPatchField, volMesh>& source,
    const labelList& faceProcAddressing,
    const labelList& cellProcAddressing,
    const labelList& boundaryProcAddressing
);

//- Map fvMatrix
template<class Type>
static void rmap
(
    fvMatrix<Type>& dest,
    const fvMatrix<Type>& source,
    const labelList& faceProcAddressing,
    const labelList& cellProcAddressing,
    const labelList& boundaryProcAddressing
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "rmapTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
