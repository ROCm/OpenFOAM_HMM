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
    Foam::upwindCPCStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef upwindCPCStencilObject_H
#define upwindCPCStencilObject_H

#include "extendedUpwindStencil.H"
#include "cellPointCellStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class upwindCPCStencilObject Declaration
\*---------------------------------------------------------------------------*/

class upwindCPCStencilObject
:
    public MeshObject<fvMesh, upwindCPCStencilObject>,
    public extendedUpwindStencil
{

public:

    TypeName("upwindCFCStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit upwindCPCStencilObject
        (
            const fvMesh& mesh,
            const scalar minOpposedness
        )
        :
            MeshObject<fvMesh, upwindCPCStencilObject>(mesh),
            extendedUpwindStencil(cellPointCellStencil(mesh), minOpposedness)
        {}


    // Destructor

        virtual ~upwindCPCStencilObject()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
