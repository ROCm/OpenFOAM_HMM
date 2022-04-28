/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam

Description
    Miscellaneous file handling for meshes.

SourceFiles
    loadOrCreateMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_loadOrCreateMesh_H
#define Foam_loadOrCreateMesh_H

#include "fvMeshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMesh;

//- Check for availability of specified mesh file (default: "faces")
boolList haveMeshFile
(
    const Time& runTime,
    const fileName& meshPath,
    const word& meshFile = "faces",
    const bool verbose = true
);


//- Remove procAddressing
void removeProcAddressing(const faMesh& mesh);

//- Remove procAddressing
void removeProcAddressing(const polyMesh& mesh);

//- Remove empty directory. Return true if removed.
bool removeEmptyDir(const fileName& path);

//- Remove empty directories from bottom up
void removeEmptyDirs(const fileName& path);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
