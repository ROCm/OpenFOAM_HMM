/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "foamVtkWriteTopoSet.H"
#include "polyMesh.H"
#include "topoSet.H"
#include "faceSet.H"
#include "cellSet.H"
#include "pointSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::vtk::writeTopoSet
(
    const polyMesh& mesh,
    const topoSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
{
    if (isA<pointSet>(set))
    {
        return vtk::writePointSet
        (
            mesh,
            dynamicCast<const pointSet&>(set),
            opts,
            file,
            parallel
        );
    }
    else if (isA<faceSet>(set))
    {
        return vtk::writeFaceSet
        (
            mesh,
            dynamicCast<const faceSet&>(set),
            opts,
            file,
            parallel
        );
    }
    else if (isA<cellSet>(set))
    {
        return vtk::writeCellSetFaces
        (
            mesh,
            dynamicCast<const cellSet&>(set),
            opts,
            file,
            parallel
        );
    }

    WarningInFunction
        << "No VTK writer for '" << set.type() << "' topoSet" << nl << endl;

    return false;
}


// ************************************************************************* //
