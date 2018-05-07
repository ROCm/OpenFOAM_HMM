/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
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

Class
    Foam::sampledTriSurfaceMeshNormal

Description
    Variant of sampledTriSurfaceMesh that samples the surface-normal component
    of a vector field.

    Returns a vector field with the value in the first component and sets
    the other two to zero.

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    (
        surface1
        {
            type    sampledTriSurfaceMeshNormal;
            surface something.obj;
            source  cells;
        }
    );
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                             | Required | Default
        type     | sampledTriSurfaceMeshNormal             | yes      |
        surface  | surface name in triSurface/             | yes      |
        source   | cells/insideCells/boundaryFaces         | yes      |
        keepIds  | pass through id numbering               | no       | false
    \endtable

SeeAlso
    Foam::sampledTriSurfaceMesh

SourceFiles
    sampledTriSurfaceMeshNormal.C
    sampledTriSurfaceMeshNormalTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledTriSurfaceMeshNormal_H
#define sampledTriSurfaceMeshNormal_H

#include "sampledTriSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sampledTriSurfaceMeshNormal Declaration
\*---------------------------------------------------------------------------*/

class sampledTriSurfaceMeshNormal
:
    public sampledTriSurfaceMesh
{
public:

    //- Runtime type information
    TypeName("sampledTriSurfaceMeshNormal");


    // Constructors

        //- Construct from components
        sampledTriSurfaceMeshNormal
        (
            const word& name,
            const polyMesh& mesh,
            const word& surfaceName,
            const samplingSource sampleSource
        );

        //- Construct from dictionary
        sampledTriSurfaceMeshNormal
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from triSurface
        sampledTriSurfaceMeshNormal
        (
            const word& name,
            const polyMesh& mesh,
            const triSurface& surface,
            const word& sampleSourceName
        );


    //- Destructor
    virtual ~sampledTriSurfaceMeshNormal() = default;


    // Member Functions

    // Sample

        //- Sample volume field onto surface faces
        virtual tmp<scalarField> sample
        (
            const interpolation<scalar>& sampler
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Sample volume field onto surface faces
        virtual tmp<vectorField> sample
        (
            const interpolation<vector>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<sphericalTensorField> sample
        (
            const interpolation<sphericalTensor>& sampler
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Sample volume field onto surface faces
        virtual tmp<symmTensorField> sample
        (
            const interpolation<symmTensor>& sampler
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Sample volume field onto surface faces
        virtual tmp<tensorField> sample
        (
            const interpolation<tensor>& sampler
        ) const
        {
            NotImplemented;
            return nullptr;
        }


    // Interpolate

        //- Interpolate volume field onto surface points
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>& interpolator
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Interpolate volume field onto surface points
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>& interpolator
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Interpolate volume field onto surface points
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>& interpolator
        ) const
        {
            NotImplemented;
            return nullptr;
        }

        //- Interpolate volume field onto surface points
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>& interpolator
        ) const
        {
            NotImplemented;
            return nullptr;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
