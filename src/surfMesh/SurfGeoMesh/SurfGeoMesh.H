/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::SurfGeoMesh

Description
    Generic mesh wrapper used by surfMesh etc.

    Similar to Foam::GeoMesh, but for meshes without an underlying
    BoundaryMesh.

\*---------------------------------------------------------------------------*/

#ifndef SurfGeoMesh_H
#define SurfGeoMesh_H

#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SurfGeoMesh Declaration
\*---------------------------------------------------------------------------*/

template<class MESH>
class SurfGeoMesh
{

protected:

    // Protected data

        //- Reference to Mesh
        const MESH& mesh_;

public:

    // Public typedefs

        typedef MESH Mesh;
        // typedef typename MESH::BoundaryMesh BoundaryMesh;


    // Constructors

        //- Construct from MESH reference
        explicit SurfGeoMesh(const MESH& mesh)
        :
            mesh_(mesh)
        {}


    // Member Functions

        //- Return the object registry
        const objectRegistry& thisDb() const
        {
            return mesh_;
        }


    // Member Operators

        //- Return reference to underlying mesh
        const MESH& operator()() const
        {
            return mesh_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
