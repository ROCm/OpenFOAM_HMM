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
    Foam::distanceSurface

Description
    A surface defined by a distance from an input searchable surface.
    Uses an iso-surface algorithm (cell, topo, point) for constructing the
    distance surface.

Usage
    Dictionary controls:
    \table
        Property | Description                              | Required | Default
        distance | distance from surface                    | yes |
        signed   | Use sign when distance is positive       | partly |
        isoMethod | Iso-algorithm (cell/topo/point)         | no  | topo
        regularise | Face simplification (enum or bool)     | no  | true
        bounds   | Limit with bounding box                  | no  |
        surfaceType | Type of surface                       | yes |
        surfaceName | Name of surface in \c triSurface/     | no  | dict name
    \endtable

Note
    For distance = 0, some special adjustments.
    - Always signed (ignoring the input value).
    - Use normal distance from surface (for better treatment of open edges).
    - When the isoSurfaceCell algorithm is used, additional checks for open
      surfaces edges are used to limit the extend of resulting distance
      surface. The resulting surface elements will not, however, contain
      partial cell coverage.

The keyword \c cell (bool value) which was use in 1906 and earlier to switch
between point/cell algorithms is now ignored (2020-12).

Changed default algorithm from cell to topo (2020-12).

SourceFiles
    distanceSurface.C

\*---------------------------------------------------------------------------*/

#ifndef distanceSurface_H
#define distanceSurface_H

#include "sampledSurface.H"
#include "searchableSurface.H"
#include "isoSurfaceCell.H"
#include "isoSurfacePoint.H"
#include "isoSurfaceTopo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class distanceSurface Declaration
\*---------------------------------------------------------------------------*/

class distanceSurface
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Surface
        const autoPtr<searchableSurface> surfPtr_;

        //- Distance value
        const scalar distance_;

        //- Signed distance
        const bool signed_;

        //- Parameters for iso-surface (algorithm, filter, mergeTol, etc)
        isoSurfaceParams isoParams_;

        //- Distance to cell centres
        autoPtr<volScalarField> cellDistancePtr_;

        //- Distance to points
        scalarField pointDistance_;

        //- Constructed iso surface (ALGO_CELL)
        autoPtr<isoSurfaceCell> isoSurfCellPtr_;

        //- Constructed iso surface (ALGO_POINT)
        autoPtr<isoSurfacePoint> isoSurfPointPtr_;

        //- Constructed iso surface (ALGO_TOPO)
        autoPtr<isoSurfaceTopo> isoSurfTopoPtr_;

public:

    //- Runtime type information
    TypeName("distanceSurface");


    // Constructors

        //- Construct from dictionary
        distanceSurface
        (
            const word& defaultSurfaceName,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from components
        distanceSurface
        (
            const polyMesh& mesh,
            const bool interpolate,
            const word& surfaceType,
            const word& surfaceName,
            const scalar distance,
            const bool useSignedDistance,
            const isoSurfaceParams& params = isoSurfaceParams()
        );


    //- Destructor
    virtual ~distanceSurface() = default;


    // Member Functions

        //- Create/recreate the distance surface
        void createGeometry();

        //- The name of the underlying searchableSurface
        const word& surfaceName() const
        {
            return surfPtr_->name();
        }

        //- The distance to the underlying searchableSurface
        scalar distance() const
        {
            return distance_;
        }


        //- The underlying surface
        const meshedSurface& surface() const
        {
            if (isoSurfCellPtr_)
            {
                return *isoSurfCellPtr_;
            }
            else if (isoSurfPointPtr_)
            {
                return *isoSurfPointPtr_;
            }
            return *isoSurfTopoPtr_;
        }


        //- The underlying surface
        meshedSurface& surface()
        {
            if (isoSurfCellPtr_)
            {
                return *isoSurfCellPtr_;
            }
            else if (isoSurfPointPtr_)
            {
                return *isoSurfPointPtr_;
            }
            return *isoSurfTopoPtr_;
        }

        //- For each face, the original cell in mesh
        const labelList& meshCells() const
        {
            if (isoSurfCellPtr_)
            {
                return isoSurfCellPtr_->meshCells();
            }
            else if (isoSurfPointPtr_)
            {
                return isoSurfPointPtr_->meshCells();
            }
            return isoSurfTopoPtr_->meshCells();
        }

        //- For each face, the original cell in mesh
        labelList& meshCells()
        {
            if (isoSurfCellPtr_)
            {
                return isoSurfCellPtr_->meshCells();
            }
            else if (isoSurfPointPtr_)
            {
                return isoSurfPointPtr_->meshCells();
            }
            return isoSurfTopoPtr_->meshCells();
        }


    // Interpolate

        //- Interpolate volume field onto surface points
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& cellValues,
            const Field<Type>& pointValues
        ) const;


    // Output

        //- Print information
        void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "distanceSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
