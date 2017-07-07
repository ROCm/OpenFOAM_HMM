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
    Foam::discreteSurface

Description
    The basis for sampling from triSurfaceMesh.
    It samples on the points/triangles of the triSurface.

    - it either samples cells or (non-coupled) boundary faces

    - 6 different modes:
        - source=cells, interpolate=false:
            finds per triangle centre the nearest cell centre and uses its value
        - source=cells, interpolate=true
            finds per triangle centre the nearest cell centre.
            Per surface point checks if this nearest cell is the one containing
            point; otherwise projects the point onto the nearest point on
            the boundary of the cell (to make sure interpolateCellPoint
            gets a valid location)

        - source=insideCells, interpolate=false:
            finds per triangle centre the cell containing it and uses its value.
            Trims triangles outside mesh.
        - source=insideCells, interpolate=true
            Per surface point interpolate cell containing it.

        - source=boundaryFaces, interpolate=false:
            finds per triangle centre the nearest point on the boundary
            (uncoupled faces only) and uses the value (or 0 if the nearest
            is on an empty boundary)
        - source=boundaryFaces, interpolate=true:
            finds per triangle centre the nearest point on the boundary
            (uncoupled faces only).
            Per surface point projects the point onto this boundary face
            (to make sure interpolateCellPoint gets a valid location)

    - since it finds a nearest per triangle each triangle is guaranteed
    to be on one processor only. So after stitching (by sampledSurfaces)
    the original surface should be complete.

SourceFiles
    discreteSurface.C
    discreteSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef discreteSurface_H
#define discreteSurface_H

#include "sampledSurface.H"
#include "triSurfaceMesh.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class treeDataFace;
class meshSearch;

/*---------------------------------------------------------------------------*\
                       Class discreteSurface Declaration
\*---------------------------------------------------------------------------*/

class discreteSurface
:
    public meshedSurface
{
public:
        //- Types of communications
        enum samplingSource
        {
            cells,
            insideCells,
            boundaryFaces
        };

private:

    //- Private typedefs for convenience
        typedef meshedSurface MeshStorage;


    // Private data

        static const Enum<samplingSource> samplingSourceNames_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Do we allow interpolation?
        const bool allowInterpolate_;

        //- Do we intend to interpolate the information?
        bool interpolate_;

        //- Surface to sample on
        const triSurfaceMesh surface_;

        //- Whether to sample internal cell values or boundary values
        const samplingSource sampleSource_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;

        //- Retain element ids/order of original surface
        bool keepIds_;

        //- List of element ids/order of the original surface,
        //  when keepIds is active.
        labelList originalIds_;

        //- Search tree for all non-coupled boundary faces
        mutable autoPtr<indexedOctree<treeDataFace>> boundaryTreePtr_;

        //- For compatibility with the meshSurf interface
        labelList zoneIds_;

        //- From local surface triangle to mesh cell/face.
        labelList sampleElements_;

        //- Local points to sample per point
        pointField samplePoints_;


    // Private Member Functions

        //- Get tree of all non-coupled boundary faces
        const indexedOctree<treeDataFace>& nonCoupledboundaryTree() const;

        bool update(const meshSearch& meshSearcher);


        //- Sample the volume field onto surface,
        //  store it (temporarily) onto the given registry
        template<class Type>
        bool sampleType
        (
            const objectRegistry& store,
            const word& fieldName
        ) const;


public:

    //- Runtime type information
    TypeName("discreteSurface");


    // Constructors

        //- Construct from components
        discreteSurface
        (
            const polyMesh& mesh,
            const word& surfaceName,
            const samplingSource sampleSource,
            const bool allowInterpolate = true
        );

        //- Construct from dictionary
        discreteSurface
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const bool allowInterpolate = true
        );

        //- Construct from triSurface
        discreteSurface
        (
            const word& name,
            const polyMesh& mesh,
            const triSurface& surface,
            const word& sampleSourceName,
            const bool allowInterpolate = true
        );


    //- Destructor
    virtual ~discreteSurface();


    // Member Functions

        //- Set new zoneIds list based on the surfZoneList information
        static void setZoneMap(const surfZoneList&, labelList& zoneIds);


      // Access

        //- Access to the underlying mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Interpolation requested for surface
        bool interpolate() const;

        //- Change interpolation request for surface.
        //  Return false if the value did not change.
        //  Use with caution.
        bool interpolate(bool b);


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


        //- Const access to per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return zoneIds_;
        }

        //- If element ids/order of the original surface are kept
        bool keepIds() const
        {
            return keepIds_;
        }

        //- List of element ids/order of the original surface,
        //  when keepIds is active.
        const labelList& originalIds() const
        {
            return originalIds_;
        }


        //- Sampling boundary values instead of cell values
        bool onBoundary() const
        {
            return sampleSource_ == boundaryFaces;
        }


        //- From local surface face to mesh cell/face.
        const labelList& sampleElements() const
        {
            return sampleElements_;
        }

        //- Local points to sample per point
        const pointField& samplePoints() const
        {
            return samplePoints_;
        }


        //- Sample the volume field onto surface,
        //  store it (temporarily) onto the given registry
        virtual bool sampleAndStore
        (
            const objectRegistry& store,
            const word& fieldName
        ) const;


        //- Sample field on surface. The tmp field is 'valid' on success.
        template<class Type>
        tmp<Field<Type>> sampleType
        (
            const word& fieldName
        ) const;


        //- Sample field on surface
        template<class Type>
        tmp<Field<Type>> sampleField
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;


        //- Interpolate field on surface
        template<class Type>
        tmp<Field<Type>> interpolateField
        (
            const interpolation<Type>&
        ) const;


        //- Write information
        virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "discreteSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
