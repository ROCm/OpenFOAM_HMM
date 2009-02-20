/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::isoSurface

Description
    A surface formed by the iso value.
    After "Regularised Marching Tetrahedra: improved iso-surface extraction",
    G.M. Treece, R.W. Prager and A.H. Gee.

    Note:
    - not possible on patches of type 'empty'. There are no values on
    'empty' patch fields so even the api would have to change
     (no volScalarField as argument). Too messy.
    - in parallel the regularisation (coarsening) always takes place
    and slightly different surfaces will be created compared to non-parallel.
    The surface will still be continuous though!
    - does tets without using cell centres/cell values. Not tested.
    - regularisation can create duplicate triangles/non-manifold surfaces.
    Current handling of those is bit ad-hoc for now and not perfect.
    - uses geometric merge with fraction of bounding box as distance.
    - triangles can be between two cell centres so constant sampling
      does not make sense.
    - does not do 2D correctly, creates non-flat iso surface.

SourceFiles
    isoSurface.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurface_H
#define isoSurface_H

#include "triSurface.H"
#include "labelPair.H"
#include "pointIndexHit.H"
#include "PackedBoolList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class isoSurface Declaration
\*---------------------------------------------------------------------------*/

class isoSurface
:
    public triSurface
{
    // Private data

        enum segmentCutType
        {
            NEARFIRST,      // intersection close to e.first()
            NEARSECOND,     //  ,,                   e.second()
            NOTNEAR         // no intersection
        };

        enum cellCutType
        {
            NOTCUT,     // not cut
            SPHERE,     // all edges to cell centre cut
            CUT         // normal cut
        };


        //- Reference to mesh
        const fvMesh& mesh_;

        //- Isosurface value
        const scalar iso_;

        //- When to merge points
        const scalar mergeDistance_;


        //- Whether face might be cut
        List<cellCutType> faceCutType_;

        //- Whether cell might be cut
        List<cellCutType> cellCutType_;

        //- Estimated number of cut cells
        label nCutCells_;

        //- For every triangle the original cell in mesh
        labelList meshCells_;

        //- For every unmerged triangle point the point in the triSurface
        labelList triPointMergeMap_;


    // Private Member Functions

        //- Get location of iso value as fraction inbetween s0,s1
        scalar isoFraction
        (
            const scalar s0,
            const scalar s1
        ) const;

        //- Check if any edge of a face is cut
        bool isEdgeOfFaceCut
        (
            const scalarField& pVals,
            const face& f,
            const bool ownLower,
            const bool neiLower
        ) const;

        //- Set faceCutType,cellCutType.
        void calcCutTypes
        (
            const volScalarField& cVals,
            const scalarField& pVals
        );

        static labelPair findCommonPoints
        (
            const labelledTri&,
            const labelledTri&
        );

        static point calcCentre(const triSurface&);

        static pointIndexHit collapseSurface
        (
            pointField& localPoints,
            DynamicList<labelledTri, 64>& localTris
        );

        void getNeighbour
        (
            const labelList& boundaryRegion,
            const volScalarField& cVals,
            const label cellI,
            const label faceI,
            scalar& nbrValue,
            point& nbrPoint
        ) const;

        //- Determine per cc whether all near cuts can be snapped to single
        //  point.
        void calcSnappedCc
        (
            const labelList& boundaryRegion,
            const volScalarField& cVals,
            const scalarField& pVals,
            DynamicList<point>& snappedPoints,
            labelList& snappedCc
        ) const;

        //- Determine per point whether all near cuts can be snapped to single
        //  point.
        void calcSnappedPoint
        (
            const PackedBoolList& isBoundaryPoint,
            const labelList& boundaryRegion,
            const volScalarField& cVals,
            const scalarField& pVals,
            DynamicList<point>& snappedPoints,
            labelList& snappedPoint
        ) const;

        //- Generate single point by interpolation or snapping
        template<class Type>
        Type generatePoint
        (
            const DynamicList<Type>& snappedPoints,

            const scalar s0,
            const Type& p0,
            const label p0Index,

            const scalar s1,
            const Type& p1,
            const label p1Index
        ) const;

        template<class Type>
        void generateTriPoints
        (
            const DynamicList<Type>& snapped,

            const scalar s0,
            const Type& p0,
            const label p0Index,

            const scalar s1,
            const Type& p1,
            const label p1Index,

            const scalar s2,
            const Type& p2,
            const label p2Index,

            const scalar s3,
            const Type& p3,
            const label p3Index,

            DynamicList<Type>& points
        ) const;

        template<class Type>
        label generateFaceTriPoints
        (
            const volScalarField& cVals,
            const scalarField& pVals,

            const GeometricField<Type, fvPatchField, volMesh>& cCoords,
            const Field<Type>& pCoords,

            const DynamicList<Type>& snappedPoints,
            const labelList& snappedCc,
            const labelList& snappedPoint,
            const label faceI,

            const scalar neiVal,
            const Type& neiPt,
            const label neiSnap,

            DynamicList<Type>& triPoints,
            DynamicList<label>& triMeshCells
        ) const;

        template<class Type>
        void generateTriPoints
        (
            const volScalarField& cVals,
            const scalarField& pVals,

            const GeometricField<Type, fvPatchField, volMesh>& cCoords,
            const Field<Type>& pCoords,

            const DynamicList<Type>& snappedPoints,
            const labelList& snappedCc,
            const labelList& snappedPoint,

            DynamicList<Type>& triPoints,
            DynamicList<label>& triMeshCells
        ) const;

        triSurface stitchTriPoints
        (
            const bool checkDuplicates,
            const List<point>& triPoints,
            labelList& triPointReverseMap,  // unmerged to merged point
            labelList& triMap               // merged to unmerged triangle
        ) const;

        //- Check single triangle for (topological) validity
        static bool validTri(const triSurface&, const label);

        //- Determine edge-face addressing
        void calcAddressing
        (
            const triSurface& surf,
            List<FixedList<label, 3> >& faceEdges,
            labelList& edgeFace0,
            labelList& edgeFace1,
            Map<labelList>& edgeFacesRest
        ) const;

        //- Determine orientation
        static void walkOrientation
        (
            const triSurface& surf,
            const List<FixedList<label, 3> >& faceEdges,
            const labelList& edgeFace0,
            const labelList& edgeFace1,
            const label seedTriI,
            labelList& flipState
        );

        //- Orient surface
        static void orientSurface
        (
            triSurface&,
            const List<FixedList<label, 3> >& faceEdges,
            const labelList& edgeFace0,
            const labelList& edgeFace1,
            const Map<labelList>& edgeFacesRest
        );

        //- Is triangle (given by 3 edges) not fully connected?
        static bool danglingTriangle
        (
            const FixedList<label, 3>& fEdges,
            const labelList& edgeFace1
        );

        //- Mark all non-fully connected triangles
        static label markDanglingTriangles
        (
            const List<FixedList<label, 3> >& faceEdges,
            const labelList& edgeFace0,
            const labelList& edgeFace1,
            const Map<labelList>& edgeFacesRest,
            boolList& keepTriangles
        );

        static triSurface subsetMesh
        (
            const triSurface& s,
            const labelList& newToOldFaces,
            labelList& oldToNewPoints,
            labelList& newToOldPoints
        );

public:

    //- Runtime type information
    TypeName("isoSurface");


    // Constructors

        //- Construct from cell values and point values. Uses boundaryField
        //  for boundary values. Requires on coupled patchfields to be set
        //  to the opposite cell value.
        isoSurface
        (
            const volScalarField& cellIsoVals,
            const scalarField& pointIsoVals,
            const scalar iso,
            const bool regularise,
            const scalar mergeTol = 1E-6    // fraction of bounding box
        );


    // Member Functions

        //- For every face original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- For every unmerged triangle point the point in the triSurface
        const labelList& triPointMergeMap() const
        {
            return triPointMergeMap_;
        }

        //- Interpolates cCoords,pCoords. Takes the original fields
        //  used to create the iso surface.
        template <class Type>
        tmp<Field<Type> > interpolate
        (
            const volScalarField& cellIsoVals,
            const scalarField& pointIsoVals,
            const GeometricField<Type, fvPatchField, volMesh>& cCoords,
            const Field<Type>& pCoords
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "isoSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
