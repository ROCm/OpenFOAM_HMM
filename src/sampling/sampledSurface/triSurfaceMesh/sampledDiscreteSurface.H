/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

        //- Sample field on faces
        template<class Type>
        tmp<Field<Type>> sampleField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        ) const;


        template<class Type>
        tmp<Field<Type>>
        interpolateField(const interpolation<Type>&) const;


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

        //- Const access to per-face zone/region information
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
        bool keepIds() const
        {
            return MeshStorage::keepIds();
        }

        //- List of element ids/order of the original surface,
        //  when keepIds is active.
        const labelList& originalIds() const
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


        //- Sample field on surface
        virtual tmp<scalarField> sample
        (
            const volScalarField&
        ) const;

        //- Sample field on surface
        virtual tmp<vectorField> sample
        (
            const volVectorField&
        ) const;

        //- Sample field on surface
        virtual tmp<sphericalTensorField> sample
        (
            const volSphericalTensorField&
        ) const;

        //- Sample field on surface
        virtual tmp<symmTensorField> sample
        (
            const volSymmTensorField&
        ) const;

        //- Sample field on surface
        virtual tmp<tensorField> sample
        (
            const volTensorField&
        ) const;


        //- Interpolate field on surface
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>&
        ) const;


        //- Interpolate field on surface
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>&
        ) const;

        //- Write information
        virtual void print(Ostream&) const;
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
