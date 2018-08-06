/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::sampledCuttingPlane

Description
    A sampledSurface defined by a plane using the iso-surface algorithm
    to \a cut the mesh.

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    (
        surface1
        {
            type        cuttingPlane;
            planeType   pointAndNormal;
            pointAndNormalDict
            {
                ...
            }
        }
    );
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                             | Required | Default
        type     | cuttingPlane                            | yes      |
        planeType | plane description (pointAndNormal etc) | yes      |
        mergeTol | tolerance for merging points            | no       | 1e-6
        regularise | point snapping                        | no       | true
        bounds   | limit with bounding box                 | no       |
        zone     | limit to cell zone (name or regex)      | no       |
        zones    | limit to cell zones (names, regexs)     | no       |
        exposedPatchName | name for zone subset            | partly   |
    \endtable

Note
    The keyword \c zones has priority over \c zone.

SeeAlso
    Foam::plane

SourceFiles
    sampledCuttingPlane.C

\*---------------------------------------------------------------------------*/

#ifndef sampledCuttingPlane_H
#define sampledCuttingPlane_H

#include "sampledSurface.H"
#include "isoSurface.H"
//#include "isoSurfaceCell.H"
#include "plane.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class sampledCuttingPlane Declaration
\*---------------------------------------------------------------------------*/

class sampledCuttingPlane
:
    public sampledSurface
{
    // Private data

        //- Plane
        const plane plane_;

        //- Optional bounding box to trim triangles against
        const boundBox bounds_;

        //- Merge tolerance
        const scalar mergeTol_;

        //- Whether to coarsen
        const bool regularise_;

        //- Whether to recalculate cell values as average of point values
        const bool average_;

        //- The zone or zones in which cutting is to occur
        wordRes zoneNames_;

        //- For zones: patch to put exposed faces into
        mutable word exposedPatchName_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;


        //- Mesh subset (optional: only used with zones)
        autoPtr<fvMeshSubset> subMeshPtr_;

        //- Distance to cell centres
        autoPtr<volScalarField> cellDistancePtr_;

        //- Distance to points
        scalarField pointDistance_;

        //- Constructed iso surface
        //autoPtr<isoSurfaceCell> isoSurfCellPtr_;
        autoPtr<isoSurface> isoSurfPtr_;


    // Private Member Functions

        //- Check and warn if bounding box does not intersect mesh or plane
        void checkBoundsIntersection
        (
            const plane& pln,
            const boundBox& meshBb
        ) const;

        //- Create iso surface
        void createGeometry();

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
    TypeName("sampledCuttingPlane");


    // Constructors

        //- Construct from dictionary
        sampledCuttingPlane
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledCuttingPlane() = default;


    // Member Functions

        //const isoSurfaceCell& surface() const
        const isoSurface& surface() const
        {
            return *isoSurfPtr_;
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

        //- Const access to per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return Foam::emptyLabelList;
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


        //- The underlying surface
        meshedSurface& surface()
        {
            return *isoSurfPtr_;
        }

        //- For each face, the original cell in mesh
        const labelList& meshCells() const
        {
            return isoSurfPtr_->meshCells();
        }

        //- For each face, the original cell in mesh
        labelList& meshCells()
        {
            return isoSurfPtr_->meshCells();
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

        //- Print information
        virtual void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledCuttingPlaneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
