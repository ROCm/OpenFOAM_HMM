/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::sampledTriSurfaceMesh

Description
    A sampledSurface from a triSurfaceMesh. It samples on the points/triangles
    of the triSurface.

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

    This is often embedded as part of a sampled surfaces function object.

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
        keepIds  | pass through id numbering               | no       | false
    \endtable

SourceFiles
    sampledTriSurfaceMesh.C
    sampledTriSurfaceMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledTriSurfaceMesh_H
#define sampledTriSurfaceMesh_H

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
                       Class sampledTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class sampledTriSurfaceMesh
:
    public sampledSurface,
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

        bool update(const meshSearch& meshSearcher);

public:

    //- Runtime type information
    TypeName("sampledTriSurfaceMesh");


    // Constructors

        //- Construct from components
        sampledTriSurfaceMesh
        (
            const word& name,
            const polyMesh& mesh,
            const word& surfaceName,
            const samplingSource sampleSource
        );

        //- Construct from dictionary
        sampledTriSurfaceMesh
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from triSurface
        sampledTriSurfaceMesh
        (
            const word& name,
            const polyMesh& mesh,
            const triSurface& surface,
            const word& sampleSourceName
        );


    //- Destructor
    virtual ~sampledTriSurfaceMesh();


    // Member Functions

        //- Set new zoneIds list based on the surfZoneList information
        static void setZoneMap
        (
            const surfZoneList& zoneLst,
            labelList& zoneIds
        );


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
            return zoneIds_;
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
            return keepIds_;
        }

        //- List of element ids/order of the original surface,
        //  when keepIds is active.
        virtual const labelList& originalIds() const
        {
            return originalIds_;
        }

        //- Sampling boundary values instead of cell values
        bool onBoundary() const
        {
            return sampleSource_ == boundaryFaces;
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
    #include "sampledTriSurfaceMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
