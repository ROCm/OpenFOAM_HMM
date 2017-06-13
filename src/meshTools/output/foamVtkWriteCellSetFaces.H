/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Write faces of cellSet to vtk polydata file.

    The data are the original cell ids

SourceFiles
    foamVtkWriteCellSetFaces.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkWriteCellSetFaces_H
#define foamVtkWriteCellSetFaces_H

#include "primitiveMesh.H"
#include "uindirectPrimitivePatch.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class primitiveMesh;
class cellSet;
class fileName;

namespace vtk
{

//- Write perimeter faces of cellset to vtk polydata file.
//  The data are the original cell ids
void writeCellSetFaces
(
    const primitiveMesh& mesh,
    const cellSet& set,
    const fileName& baseName,
    const vtk::outputOptions outOpts
);

} // End namespace vtk
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
