/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::sampledDiscreteSurface

Description
    A sampledSurface from a triSurfaceMesh.
    It samples on the points/triangles of a triSurfaceMesh.

See Also
    discreteSurface, sampledSurface

SourceFiles
    sampledDiscreteSurface.C
    sampledDiscreteSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledDiscreteSurface_H
#define sampledDiscreteSurface_H

#include "sampledSurface.H"
#include "discreteSurface.H"
#include "triSurfaceMesh.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sampledDiscreteSurface;

/*---------------------------------------------------------------------------*\
                   Class sampledDiscreteSurface Declaration
\*---------------------------------------------------------------------------*/

class sampledDiscreteSurface
:
    public sampledSurface,
    public discreteSurface
{
    //- Private typedefs for convenience
    typedef discreteSurface MeshStorage;
    typedef discreteSurface SurfaceSource;


    // Private Member Functions

        //- Sample volume field onto surface faces
        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler
        ) const;

        //- Interpolate volume field onto surface points
        template<class Type>
        tmp<Field<Type>> sampleOnPoints
        (
            const interpolation<Type>& interpolator
        ) const;


public:

    //- Runtime type information
    TypeName("sampledDiscreteSurface");


    // Constructors

        //- Construct from components
        sampledDiscreteSurface
        (
            const word& name,
            const polyMesh& mesh,
            const word& surfaceName,
            const discreteSurface::samplingSource sampleSource
        );

        //- Construct from dictionary
        sampledDiscreteSurface
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from triSurface
        sampledDiscreteSurface
        (
            const word& name,
            const polyMesh& mesh,
            const triSurface& surface,
            const word& sampleSourceName
        );


    //- Destructor
    virtual ~sampledDiscreteSurface();


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
        bool update(const treeBoundBox&);

        //- Points of surface
        virtual const pointField& points() const
        {
            return MeshStorage::points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return MeshStorage::surfFaces();
        }

        //- Per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return MeshStorage::zoneIds();
        }

        //- Face area vectors
        virtual const vectorField& Sf() const
        {
            return MeshStorage::Sf();
        }

        //- Face area magnitudes
        virtual const scalarField& magSf() const
        {
            return MeshStorage::magSf();
        }

        //- Face centres
        virtual const vectorField& Cf() const
        {
            return MeshStorage::Cf();
        }

        //- If element ids/order of the original surface are kept
        virtual bool hasFaceIds() const
        {
            return MeshStorage::hasFaceIds();
        }

        //- List of element ids/order of the original surface,
        //  when keepIds is active.
        virtual const labelList& originalIds() const
        {
            return MeshStorage::originalIds();
        }

        //- Sample the volume field onto surface,
        //  store it (temporarily) onto the given registry
        virtual bool sampleAndStore
        (
            const objectRegistry& store,
            const word& fieldName
        ) const;


    // Sample

        //- Sample volume field onto surface faces
        virtual tmp<scalarField> sample
        (
            const interpolation<scalar>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<vectorField> sample
        (
            const interpolation<vector>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<sphericalTensorField> sample
        (
            const interpolation<sphericalTensor>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<symmTensorField> sample
        (
            const interpolation<symmTensor>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<tensorField> sample
        (
            const interpolation<tensor>& sampler
        ) const;


    // Interpolate

        //- Interpolate volume field onto surface points
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>& interpolator
        ) const;


    // Output

        //- Write
        virtual void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "sampledDiscreteSurfaceTemplates.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
