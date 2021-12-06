/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

Class
    Foam::glTF::scene

Description
    Main class to assemble glTF components into a scene

Note
    Implements the glTF v2 specification

SourceFiles
    foamGltfScene.C
    foamGltfSceneTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foam_gltf_scene_H
#define foam_gltf_scene_H

#include "foamGltfList.H"
#include "foamGltfObject.H"
#include "foamGltfMesh.H"
#include "foamGltfBufferView.H"
#include "foamGltfAccessor.H"
#include "foamGltfAnimation.H"
#include "scalarField.H"
#include "vectorField.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace glTF
{

/*---------------------------------------------------------------------------*\
                         Class glTF::scene Declaration
\*---------------------------------------------------------------------------*/

class scene
{
    // Private Data

        //- List of binary objects
        glTF::List<object> objects_;

        //- List of meshes
        glTF::List<mesh> meshes_;

        //- List of buffer views
        glTF::List<bufferView> bufferViews_;

        //- List of accessors
        glTF::List<accessor> accessors_;

        //- List of animations
        glTF::List<animation> animations_;

        //- Accumulative size in bytes
        label bytes_;


public:

    // Constructors

        //- Default construct
        scene();


    // Public Member Functions

        //- Returns accessor index
        template<class Type>
        label addField
        (
            const Type& fld,
            const word& name,
            const label target = -1
        );

        //- Returns index of last mesh
        template<class Type>
        label addMesh(const Type& fld, const word& name);

        //- Returns accessor index
        template<class Type>
        label addFieldToMesh
        (
            const Type& fld,
            const word& name,
            const label meshi
        );

        //- Returns accessor index
        label addColourToMesh
        (
            const vectorField& fld,
            const word& name,
            const label meshi,
            const scalarField& alpha = scalarField()
        );

        //- Returns index of last animation
        label createAnimation(const word& name);

        //- Add to existing animation
        void addToAnimation
        (
            const label animationi,
            const label inputId,
            const label outputId,
            const label meshId,
            const string& interpolation = "LINEAR"
        );

        //- Write to stream (JSON and binary data)
        void write(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace glTF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamGltfSceneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
