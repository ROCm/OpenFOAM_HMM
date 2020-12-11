/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::isoSurfaceTopo

Description
    Marching tet iso surface algorithm with optional filtering to keep only
    points originating from mesh edges.

SourceFiles
    isoSurfaceTopo.C
    isoSurfaceTopoTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurfaceTopo_H
#define isoSurfaceTopo_H

#include "bitSet.H"
#include "edgeList.H"
#include "isoSurfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class isoSurfaceTopo Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceTopo
:
    public isoSurfaceBase
{
    // Private Data

        //- Corrected version of tetBasePtIs
        labelList tetBasePtIs_;

        //- Per point: originating mesh vertex/cc. See encoding above
        edgeList pointToVerts_;

        //- For every point the originating face in mesh
        labelList pointToFace_;

        //- The cell cut type
        List<cutType> cellCutType_;


    // Private Member Functions

        scalar minTetQ
        (
            const label facei,
            const label faceBasePtI
        ) const;

        void fixTetBasePtIs();

        //- Generate single point on edge
        label generatePoint
        (
            const label facei,
            const bool edgeIsDiag,
            const edge& vertices,

            DynamicList<edge>& pointToVerts,
            DynamicList<label>& pointToFace,
            DynamicList<bool>& pointFromDiag,
            EdgeMap<label>& vertsToPoint
        ) const;

        //- Generate triangles from tet
        void generateTriPoints
        (
            const label facei,
            const int tetCutIndex,  //!< Encoded tet cuts. getTetCutIndex()
            const tetCell& tetLabels,
            const FixedList<bool, 6>& edgeIsDiag,

            DynamicList<edge>& pointToVerts,
            DynamicList<label>& pointToFace,
            DynamicList<bool>& pointFromDiag,

            EdgeMap<label>& vertsToPoint,
            DynamicList<label>& verts
        ) const;

        //- Generate triangles from cell
        void generateTriPoints
        (
            const label celli,
            const bool isTet,

            DynamicList<edge>& pointToVerts,
            DynamicList<label>& pointToFace,
            DynamicList<bool>& pointFromDiag,

            EdgeMap<label>& vertsToPoint,
            DynamicList<label>& verts,
            DynamicList<label>& faceLabels
        ) const;


        // Simplification

            static void triangulateOutside
            (
                const bool filterDiag,
                const primitivePatch& pp,

                const boolUList& pointFromDiag,
                const labelUList& pointToFace,
                const label cellID,

                DynamicList<face>& compactFaces,
                DynamicList<label>& compactCellIDs
            );

            static void removeInsidePoints
            (
                Mesh& s,                            // Modify inplace
                const bool filterDiag,

                // Inputs
                const boolUList& pointFromDiag,
                const labelUList& pointToFace,
                const labelUList& start,            // Per cell:starting tri

                // Outputs
                DynamicList<label>& pointCompactMap, // Per point the original
                DynamicList<label>& compactCellIDs   // Per face the cellID
            );


protected:

    // Editing

        //- Subset the surface using the selected faces.
        //
        //  \param[in] include the faces to select
        void inplaceSubsetMesh(const bitSet& include);


    // Sampling

        //- Interpolates cCoords,pCoords.
        template<class Type>
        tmp<Field<Type>> interpolateTemplate
        (
            const Field<Type>& cCoords,
            const Field<Type>& pCoords
        ) const;

public:

    //- Runtime type information
    TypeName("isoSurfaceTopo");


    // Constructors

        //- Construct from cell and point values
        //
        //  \param ignoreCells cells to ignore in the cellValues
        //
        //  Control parameters include
        //  - bounds optional bounding box for trimming
        //  - mergeTol fraction of mesh bounding box for merging points
        isoSurfaceTopo
        (
            const polyMesh& mesh,
            const scalarField& cellValues,
            const scalarField& pointValues,
            const scalar iso,
            const isoSurfaceParams& params = isoSurfaceParams(),
            const bitSet& ignoreCells = bitSet()
        );


    //- Destructor
    virtual ~isoSurfaceTopo() = default;


    // Member Functions

        //- For every point originating face (pyramid) in mesh
        const labelList& pointToFace() const
        {
            return pointToFace_;
        }

        //- Per point: originating mesh vertex/cc. See encoding above
        const edgeList& pointToVerts() const
        {
            return pointToVerts_;
        }


    // Sampling

        declareIsoSurfaceInterpolateMethod(scalar);
        declareIsoSurfaceInterpolateMethod(vector);
        declareIsoSurfaceInterpolateMethod(sphericalTensor);
        declareIsoSurfaceInterpolateMethod(symmTensor);
        declareIsoSurfaceInterpolateMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "isoSurfaceTopoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
