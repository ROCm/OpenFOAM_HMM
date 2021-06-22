/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::fvMeshSubset

Description
    Given the original mesh and the list of selected cells, it creates the
    mesh consisting only of the desired cells, with the mapping list for
    points, faces, and cells.

    Puts all exposed internal faces into either
    - a user supplied patch
    - a newly created patch "oldInternalFaces"

    - setCellSubset does coupled patch subsetting as well. If it detects
      a face on a coupled patch 'losing' its neighbour it will move the
      face into the oldInternalFaces patch.

    - if a user supplied patch is used it is up to the destination
      patchField to handle exposed internal faces (mapping from face -1).
      If not provided the default is to assign the internalField. All the
      basic patch field types (e.g. fixedValue) will give a warning and
      preferably derived patch field types should be used that know how to
      handle exposed faces (e.g. use uniformFixedValue instead of fixedValue)

SourceFiles
    fvMeshSubset.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshSubset_H
#define fvMeshSubset_H

#include "fvMesh.H"
#include "pointMesh.H"
#include "GeometricField.H"
#include "bitSet.H"
#include "HashSet.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fvMeshSubset Declaration
\*---------------------------------------------------------------------------*/

class fvMeshSubset
{
    // Private Data

        //- Mesh to subset from
        const fvMesh& baseMesh_;

        //- Demand-driven subset mesh pointer
        autoPtr<fvMesh> fvMeshSubsetPtr_;

        //- Optional face mapping array with flip encoded (-1/+1)
        mutable autoPtr<labelList> faceFlipMapPtr_;

        //- Point mapping array
        labelList pointMap_;

        //- Face mapping array
        labelList faceMap_;

        //- Cell mapping array
        labelList cellMap_;

        //- Patch mapping array
        labelList patchMap_;


    // Private Member Functions

        //- Check if subset has been performed
        bool checkCellSubset() const;

        //- Calculate face flip map
        void calcFaceFlipMap() const;

        //- Adapt nCellsUsingFace for coupled faces becoming 'uncoupled'.
        void doCoupledPatches
        (
            const bool syncPar,
            labelList& nCellsUsingFace
        ) const;

        //- Forwarding to Foam::removeCells
        void removeCellsImpl
        (
            const bitSet& cellsToRemove,
            const labelList& exposedFaces,
            const labelList& patchIDs,
            const bool syncCouples
        );


        //- Subset of subset
        static labelList subsetSubset
        (
            const label nElems,
            const labelUList& selectedElements,
            const labelUList& subsetMap
        );

        //- Create zones for submesh
        void subsetZones();

        //- Helper: extract cells-to-remove from cells-to-keep
        bitSet getCellsToRemove
        (
            const bitSet& selectedCells
        ) const;

        //- Helper: extract cells-to-remove from cells-to-keep
        bitSet getCellsToRemove
        (
            const label regioni,
            const labelUList& regions
        ) const;


        //- No copy construct
        fvMeshSubset(const fvMeshSubset&) = delete;

        //- No copy assignment
        void operator=(const fvMeshSubset&) = delete;


public:

    // Static Data Members

        //- Name for exposed internal faces (default: oldInternalFaces)
        static word exposedPatchName;


    // Constructors

        //- Construct given a mesh to subset
        explicit fvMeshSubset(const fvMesh& baseMesh);

        //- Construct for a cell-subset of the given mesh
        //  See setCellSubset() for more details.
        fvMeshSubset
        (
            const fvMesh& baseMesh,
            const bitSet& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Construct for a cell-subset of the given mesh
        //  See setCellSubset() for more details.
        fvMeshSubset
        (
            const fvMesh& baseMesh,
            const labelUList& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Construct for a cell-subset of the given mesh
        //  See setCellSubset() for more details.
        fvMeshSubset
        (
            const fvMesh& baseMesh,
            const labelHashSet& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Construct for a cell-subset of the given mesh
        //  See setCellSubset() for more details.
        fvMeshSubset
        (
            const fvMesh& baseMesh,
            const label regioni,
            const labelUList& regions,
            const label patchID = -1,
            const bool syncPar = true
        );


    // Member Functions

    // Access

        //- Original mesh
        inline const fvMesh& baseMesh() const noexcept;

        //- Return baseMesh or subMesh, depending on the current state.
        inline const fvMesh& mesh() const noexcept;

        //- Have subMesh?
        inline bool hasSubMesh() const noexcept;

        //- Return reference to subset mesh
        inline const fvMesh& subMesh() const;

        //- Return reference to subset mesh
        inline fvMesh& subMesh();

        //- Return point map
        inline const labelList& pointMap() const;

        //- Return face map
        inline const labelList& faceMap() const;

        //- Return face map with sign to encode flipped faces
        inline const labelList& faceFlipMap() const;

        //- Return cell map
        inline const labelList& cellMap() const;

        //- Return patch map
        inline const labelList& patchMap() const;


    // Edit

        //- Reset maps and subsetting
        void clear();


        //- Define cell subset based on the selectedCells.
        //  Create "oldInternalFaces" patch for exposed
        //  internal faces (patchID==-1) or use supplied patch.
        //  Handles coupled patches if necessary by making coupled patch
        //  face part of patchID (so uncoupled)
        void setCellSubset
        (
            const bitSet& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Define cell subset, using the specified cells
        //- to define the selection
        void setCellSubset
        (
            const labelUList& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Define cell subset, using the specified cells
        //- labelHashSet to define the selection
        void setCellSubset
        (
            const labelHashSet& selectedCells,
            const label patchID = -1,
            const bool syncPar = true
        );

        //- Define cell subset, using the cells for which region == regioni.
        void setCellSubset
        (
            const label regioni,
            const labelUList& regions,
            const label patchID = -1,
            const bool syncCouples = true
        );


    // Two-step subsetting

        //- Get labels of exposed faces.
        //  These are
        //  - internal faces that become boundary faces
        //  - coupled faces that become uncoupled (since one of the
        //    sides gets deleted)
        labelList getExposedFaces
        (
            const bitSet& selectedCells,
            const bool syncCouples = true
        ) const;

        //- Get labels of exposed faces.
        //  These are
        //  - internal faces that become boundary faces
        //  - coupled faces that become uncoupled (since one of the
        //    sides gets deleted)
        labelList getExposedFaces
        (
            const label regioni,
            const labelUList& regions,
            const bool syncCouples = true
        ) const;

        //- For every exposed face (from above getExposedFaces)
        //  used supplied (existing!) patch
        void setCellSubset
        (
            const bitSet& selectedCells,
            const labelList& exposedFaces,
            const labelList& patchIDs,
            const bool syncCouples = true
        );

        //- For every exposed face (from above getExposedFaces)
        //  used supplied (existing!) patch
        void setCellSubset
        (
            const label regioni,
            const labelList& regions,
            const labelList& exposedFaces,
            const labelList& patchIDs,
            const bool syncCouples = true
        );


    // Field Mapping

        //- Map volume field. Optionally allow unmapped faces not to produce
        //  a warning
        template<class Type>
        static tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const fvMesh& sMesh,
            const labelUList& patchMap,
            const labelUList& cellMap,
            const labelUList& faceMap,
            const bool allowUnmapped = false
        );

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const bool allowUnmapped = false
        ) const;

        //- Map surface field. Optionally negates value if flipping
        //  a face (from exposing an internal face)
        template<class Type>
        static tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&,
            const fvMesh& sMesh,
            const labelUList& patchMap,
            const labelUList& cellMap,
            const labelUList& faceMap
        );

        //- Map surface field. Optionally allow unmapped faces not to produce
        //  a warning (not currently used)
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&,
            const bool allowUnmapped = false
        ) const;

        //- Map point field
        template<class Type>
        static tmp<GeometricField<Type, pointPatchField, pointMesh>>
        interpolate
        (
            const GeometricField<Type, pointPatchField, pointMesh>&,
            const pointMesh& sMesh,
            const labelUList& patchMap,
            const labelUList& pointMap
        );

        //- Map point field. Optionally allow unmapped points not to produce
        //  a warning (not currently used)
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>>
        interpolate
        (
            const GeometricField<Type, pointPatchField, pointMesh>&,
            const bool allowUnmapped = false
        ) const;

        //- Map dimensioned field
        template<class Type>
        static tmp<DimensionedField<Type, volMesh>>
        interpolate
        (
            const DimensionedField<Type, volMesh>&,
            const fvMesh& sMesh,
            const labelUList& cellMap
        );

        //- Map Dimensioned. Optional unmapped argument not used
        template<class Type>
        tmp<DimensionedField<Type, volMesh>>
        interpolate
        (
            const DimensionedField<Type, volMesh>&,
            const bool allowUnmapped = false
        ) const;


    // Compatibility

        //- Deprecated(2018-07) old method name and old parameter ordering.
        //  \deprecated(2018-07) - use setCellSubset() method
        inline void setLargeCellSubset
        (
            const labelUList& region,
            const label currentRegion,
            const label patchID = -1,
            const bool syncCouples = true
        )
        {
            Info<< "WARNING: using highly deprecated method: "
                << "fvMeshSubset::setLargeCellSubset()" << nl;

            setCellSubset
            (
                currentRegion,
                region,
                patchID,
                syncCouples
            );
        }


        //- Deprecated(2018-07) old method name
        //  \deprecated(2018-07) - use setCellSubset() method
        inline void setLargeCellSubset
        (
            const labelHashSet& globalCellMap,
            const label patchID = -1,
            const bool syncPar = true
        )
        {
            Info<< "WARNING: using highly deprecated method: "
                << "fvMeshSubset::setLargeCellSubset()" << nl;

            setCellSubset(globalCellMap, patchID, syncPar);
        }


        //- Deprecated(2018-07) method
        //  For every exposed face (from getExposedFaces) use supplied
        //  (existing!) patch ids
        //  \deprecated(2018-07) - use setCellSubset() method
        inline void setLargeCellSubset
        (
            const labelList& regions,
            const label regioni,
            const labelList& exposedFaces,
            const labelList& patchIDs,
            const bool syncCouples = true
        )
        {
            Info<< "WARNING: using highly deprecated method: "
                << "fvMeshSubset::setLargeCellSubset()" << nl;

            setCellSubset
            (
                regioni,
                regions,
                exposedFaces,
                patchIDs,
                syncCouples
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMeshSubsetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshSubsetInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
