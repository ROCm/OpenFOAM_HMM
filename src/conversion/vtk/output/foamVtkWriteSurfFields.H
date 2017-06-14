/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

InNamespace
    Foam::vtk

Description
    Write surface fields as vectors

SourceFiles
    foamVtkWriteSurfFields.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkWriteSurfFields_H
#define foamVtkWriteSurfFields_H

#include "fvMesh.H"
#include "surfaceMesh.H"
#include "surfaceFieldsFwd.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

//- Write surface vector fields
void writeSurfFields
(
    const fvMesh& mesh,
    const fileName& baseName,
    const vtk::outputOptions outOpts,
    const UPtrList<const surfaceVectorField>& surfVectorFields
);

} // End namespace vtk
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
