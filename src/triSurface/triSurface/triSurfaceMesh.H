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

Class
    Foam::triSurfaceMesh

Description
    A triSurface with objectRegistry so it can hold fields.

SourceFiles
    triSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMesh_H
#define triSurfaceMesh_H

#include "triSurface.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class triSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceMesh
:
    public objectRegistry,
    public triSurface
{

public:

    // Constructors

        //- Construct from triangles, patches, points.
        triSurfaceMesh
        (
            const IOobject& io,
            const List<labelledTri>& faces,
            const geometricSurfacePatchList& patches,
            const pointField& points
        )
        :
            objectRegistry(io),
            triSurface(faces, patches, points)
        {}

        //- Construct from file name (uses extension to determine type)
        triSurfaceMesh(const IOobject& io, const fileName& fName)
        :
            objectRegistry(io),
            triSurface(fName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
