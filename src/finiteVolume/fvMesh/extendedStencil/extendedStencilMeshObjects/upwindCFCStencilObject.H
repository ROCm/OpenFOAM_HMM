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
    Foam::upwindCFCStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef upwindCFCStencilObject_H
#define upwindCFCStencilObject_H

#include "extendedUpwindStencil.H"
#include "cellFaceCellStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class upwindCFCStencilObject Declaration
\*---------------------------------------------------------------------------*/

class upwindCFCStencilObject
:
    public MeshObject<fvMesh, upwindCFCStencilObject>,
    public extendedUpwindStencil
{

public:

    TypeName("upwindCFCStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit upwindCFCStencilObject
        (
            const fvMesh& mesh,
            const scalar minOpposedness
        )
        :
            MeshObject<fvMesh, upwindCFCStencilObject>(mesh),
            extendedUpwindStencil(cellFaceCellStencil(mesh), minOpposedness)
        {}


    // Destructor

        virtual ~upwindCFCStencilObject()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
