/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::sampledIsoSurfaceTopo

Description
    A sampledSurface defined by a surface of iso value.
    To be used in sampleSurfaces / functionObjects. Recalculates iso surface
    only if time changes.

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    (
        surface1
        {
            type    isoSurfaceTopo;
            isoField        p;
            isoValue        0.0;
        }
    );
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                             | Required | Default
        type     | isoSurfaceTopo                          | yes      |
        isoField | field name for obtaining iso-surface    | yes      |
        isoValue | value of iso-surface                    | yes      |
        regularise | filter faces                          | no       | true
        triangulate | triangulate faces (if regularise)    | no       | false
    \endtable

Note
    Does not currently support cell zones.

SourceFiles
    sampledIsoSurfaceTopo.C

\*---------------------------------------------------------------------------*/

#ifndef sampledIsoSurfaceTopo_H
#define sampledIsoSurfaceTopo_H

#include "sampledSurface.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class sampledIsoSurfaceTopo Declaration
\*---------------------------------------------------------------------------*/

class sampledIsoSurfaceTopo
:
    public sampledSurface,
    public MeshedSurface<face>
{
    // Private typedefs for convenience
    typedef MeshedSurface<face> MeshStorage;

    // Private data

        //- Field to get isoSurface of
        const word isoField_;

        //- Iso value
        const scalar isoVal_;

        //- Whether to coarse
        const bool regularise_;

        //- Whether to triangulate
        const bool triangulate_;

    // Recreated for every isoSurface

        //- Time at last call, also track it surface needs an update
        mutable label prevTimeIndex_;

        //- For every triangle/face the original cell in mesh
        mutable labelList meshCells_;


    // Private Member Functions

        //- Create iso surface (if time has changed)
        //  Do nothing (and return false) if no update was needed
        bool updateGeometry() const;

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
    TypeName("sampledIsoSurfaceTopo");


    // Constructors

        //- Construct from dictionary
        sampledIsoSurfaceTopo
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledIsoSurfaceTopo();


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


        //- Points of surface
        virtual const pointField& points() const
        {
            return MeshStorage::points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return *this;
        }

        //- Per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return labelList::null();
        }

        //- Face area magnitudes
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

        //- Write
        virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledIsoSurfaceTopoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
