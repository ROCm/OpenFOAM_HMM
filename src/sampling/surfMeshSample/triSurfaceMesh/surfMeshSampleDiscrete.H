/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::surfMeshSampleDiscrete

Description
    Sampling surfFields onto a surfMesh based on a triSurfaceMesh.

    This is often embedded as part of a surfMeshSamplers function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    (
        surface1
        {
            type    sampledTriSurfaceMesh;
            surface something.obj;
            source  cells;
        }
    );
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                             | Required | Default
        type     | sampledTriSurfaceMesh                   | yes      |
        surface  | surface name in triSurface/             | yes      |
        source   | cells/insideCells/boundaryFaces         | yes      |
    \endtable

See Also
    discreteSurface, surfMeshSample

SourceFiles
    surfMeshSampleDiscrete.C

\*---------------------------------------------------------------------------*/

#ifndef surfMeshSampleDiscrete_H
#define surfMeshSampleDiscrete_H

#include "surfMeshSample.H"
#include "discreteSurface.H"
#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class surfMeshSampleDiscrete;

/*---------------------------------------------------------------------------*\
                       Class surfMeshSampleDiscrete Declaration
\*---------------------------------------------------------------------------*/

class surfMeshSampleDiscrete
:
    public surfMeshSample,
    private discreteSurface
{
    // Private typedefs for convenience
    typedef discreteSurface SurfaceSource;


    // Private Member Functions

        //- Transfer mesh content from SurfaceSource to surfMesh
        void transferContent();

        //- Sample field (from volume field) onto surface faces
        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler
        ) const;

        //- Sample field on surface.
        template<class Type>
        bool sampleType
        (
            const word& fieldName,
            const word& sampleScheme
        ) const;


public:

    //- Runtime type information
    TypeName("surfMeshSampleDiscrete");


    // Constructors

        //- Construct from components
        surfMeshSampleDiscrete
        (
            const word& name,
            const polyMesh& mesh,
            const word& surfaceName,
            const discreteSurface::samplingSource sampleSource
        );

        //- Construct from dictionary
        surfMeshSampleDiscrete
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from triSurface
        surfMeshSampleDiscrete
        (
            const word& name,
            const polyMesh& mesh,
            const triSurface& surface,
            const word& sampleSourceName
        );


    //- Destructor
    virtual ~surfMeshSampleDiscrete() = default;


    // Member Functions

        //- Does the surface need an update?
        virtual bool needsUpdate() const;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire();

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was needed
        virtual bool update();

        //- Update the surface using a bound box to limit the searching.
        //  For direct use, i.e. not through sample.
        //  Do nothing (and return false) if no update was needed
        bool update(const treeBoundBox& bb);

        //- Sample the volume field onto surface
        virtual bool sample
        (
            const word& fieldName,
            const word& sampleScheme = "cell"
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfMeshSampleDiscreteTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
