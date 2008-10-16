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
    Foam::centredCECStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCECStencilObject_H
#define centredCECStencilObject_H

#include "extendedCentredStencil.H"
#include "cellEdgeCellStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class centredCECStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCECStencilObject
:
    public MeshObject<fvMesh, centredCECStencilObject>,
    public extendedCentredStencil
{

public:

    TypeName("centredCECStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit centredCECStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject<fvMesh, centredCECStencilObject>(mesh),
            extendedCentredStencil(cellEdgeCellStencil(mesh))
        {}


    // Destructor

        virtual ~centredCECStencilObject()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
