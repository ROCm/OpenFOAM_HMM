/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::sampledIsoSurface

Description
    A sampledSurface defined by a surface of iso value using a
    \em point algorithm (always triangulated!).
    It only recalculates the iso-surface if time changes.
    To be used in sampleSurfaces / functionObjects.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type     isoSurfacePoint;
            isoField T;
            isoValue 373;
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                              | Required | Default
        type     | isoSurfacePoint / isoSurface             | yes |
        isoField | field name for obtaining iso-surface     | yes |
        isoValue | value of iso-surface                     | yes |
        average  | cell values from averaged point values   | no  | false
        bounds   | limit with bounding box                  | no  |
        zone     | limit to cell zone (name or regex)       | no  |
        zones    | limit to cell zones (names, regexs)      | no  |
        exposedPatchName | name for zone subset             | partly |
        regularise | point snapping (bool or enum)          | no  | true
        mergeTol | tolerance for merging points             | no  | 1e-6
    \endtable

SourceFiles
    sampledIsoSurface.C
    sampledIsoSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledIsoSurface_H
#define sampledIsoSurface_H

#include "isoSurfacePoint.H"
#include "sampledSurface.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class sampledIsoSurface Declaration
\*---------------------------------------------------------------------------*/

class sampledIsoSurface
:
    public sampledSurface
{
    // Private Data

        //- Field to get isoSurface of
        const word isoField_;

        //- Iso value
        const scalar isoVal_;

        //- Parameters (filtering etc) for iso-surface
        isoSurfaceParams isoParams_;

        //- Whether to recalculate cell values as average of point values
        bool average_;

        //- The zone or zones for the iso-surface
        wordRes zoneNames_;

        //- For zones: patch to put exposed faces into
        mutable word exposedPatchName_;

        mutable autoPtr<isoSurfacePoint> isoSurfacePtr_;


    // Recreated for every isoSurface

        //- Time at last call, also track if surface needs an update
        mutable label prevTimeIndex_;

        //- Cached volfield
        mutable autoPtr<volScalarField> storedVolFieldPtr_;
        mutable const volScalarField* volFieldPtr_;

        //- Cached pointfield
        mutable const pointScalarField* pointFieldPtr_;


    // And on subsetted mesh

        //- Cached submesh
        mutable autoPtr<fvMeshSubset> subMeshPtr_;

        //- Cached volfield
        mutable autoPtr<volScalarField> storedVolSubFieldPtr_;
        mutable const volScalarField* volSubFieldPtr_;

        //- Cached pointfield
        mutable const pointScalarField* pointSubFieldPtr_;


    // Private Member Functions

        //- Get fields needed to recreate iso surface.
        void getIsoFields() const;

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
    TypeName("samplediIsoSurfacePoint");


    // Constructors

        //- Construct from dictionary
        sampledIsoSurface
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledIsoSurface();


    // Member Functions

        const isoSurfacePoint& surface() const
        {
            return *isoSurfacePtr_;
        }

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
            return surface().points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return surface().surfFaces();
        }

        //- Per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return labelList::null();
        }

        //- Face area magnitudes
        virtual const vectorField& Sf() const
        {
            return surface().Sf();
        }

        //- Face area magnitudes
        virtual const scalarField& magSf() const
        {
            return surface().magSf();
        }

        //- Face centres
        virtual const vectorField& Cf() const
        {
            return surface().Cf();
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


    // Output

        //- Write
        virtual void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledIsoSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
