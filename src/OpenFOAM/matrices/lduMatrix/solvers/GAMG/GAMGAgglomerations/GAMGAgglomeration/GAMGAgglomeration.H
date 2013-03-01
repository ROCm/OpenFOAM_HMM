/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::GAMGAgglomeration

Description
    Geometric agglomerated algebraic multigrid agglomeration class.

SourceFiles
    GAMGAgglomeration.C
    GAMGAgglomerationTemplates.C
    GAMGAgglomerate.C
    GAMGAgglomerateLduAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef GAMGAgglomeration_H
#define GAMGAgglomeration_H

#include "MeshObject.H"
#include "lduPrimitiveMesh.H"
#include "lduInterfacePtrsList.H"
#include "primitiveFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduMesh;
class lduMatrix;
class mapDistribute;

/*---------------------------------------------------------------------------*\
                    Class GAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class GAMGAgglomeration
:
    public MeshObject<lduMesh, GAMGAgglomeration>
{
protected:

    // Protected data

        //- Max number of levels
        label maxLevels_;

        //- Number of cells in coarsest level
        label nCellsInCoarsestLevel_;

        //- The number of cells in each level
        labelList nCells_;

        //- Cell restriction addressing array.
        //  Maps from the finer to the coarser level.
        PtrList<labelField> restrictAddressing_;

        //- For cross-processor operation
        autoPtr<mapDistribute> restrictDistributeMapPtr_;

        //- Face restriction addressing array.
        //  Maps from the finer to the coarser level.
        //  Positive indices map the finer faces which form part of the boundary
        //  of the coarser cells to the corresponding coarser cell face.
        //  Negative indices map the finer faces which are internal to the
        //  coarser cells to minus the corresponding coarser cell index minus 1.
        PtrList<labelList> faceRestrictAddressing_;

        //- For cross-processor operation
        autoPtr<mapDistribute> faceRestrictDistributeMapPtr_;

        //- Hierarchy of mesh addressing
        PtrList<lduPrimitiveMesh> meshLevels_;

        //- Hierarchy interfaces.
        //  Warning: Needs to be deleted explicitly.
        PtrList<lduInterfacePtrsList> interfaceLevels_;


        //- Processor agglomeration

            //- Combined coarsest meshes for all processors
            mutable autoPtr<lduPrimitiveMesh> allMeshPtr_;

            //- Per remote processor the local equivalent. TBD: hack. Can be
            //  deleted once lduMatrices get received and merged into single
            //  one without constructing local ones
            mutable PtrList<lduMesh> otherMeshes_;

            //- Mapping from processor to coarsestAllMesh cells
            mutable labelList cellOffsets_;

            //- Mapping from processor to coarsestAllMesh face
            mutable labelListList faceMap_;

            //- Mapping from processor to coarsestAllMesh boundary
            mutable labelListList boundaryMap_;

            //- Mapping from processor to coarsestAllMesh boundary face
            mutable labelListListList boundaryFaceMap_;


        //- Assemble coarse mesh addressing
        void agglomerateLduAddressing(const label fineLevelIndex);

        //- Shrink the number of levels to that specified
        void compactLevels(const label nCreatedLevels);

        //- Check the need for further agglomeration
        bool continueAgglomerating(const label nCoarseCells) const;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        GAMGAgglomeration(const GAMGAgglomeration&);

        //- Disallow default bitwise assignment
        void operator=(const GAMGAgglomeration&);


public:

    //- Runtime type information
    TypeName("GAMGAgglomeration");


    // Declare run-time constructor selection tables

        //- Runtime selection table for pure geometric agglomerators
        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGAgglomeration,
            lduMesh,
            (
                const lduMesh& mesh,
                const dictionary& controlDict
            ),
            (
                mesh,
                controlDict
            )
        );

        //- Runtime selection table for matrix or mixed geometric/matrix
        //  agglomerators
        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGAgglomeration,
            lduMatrix,
            (
                const lduMatrix& matrix,
                const dictionary& controlDict
            ),
            (
                matrix,
                controlDict
            )
        );


    // Constructors

        //- Construct given mesh and controls
        GAMGAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );


    // Selectors

        //- Return the selected geometric agglomerator
        static const GAMGAgglomeration& New
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );

        //- Return the selected matrix agglomerator
        static const GAMGAgglomeration& New
        (
            const lduMatrix& matrix,
            const dictionary& controlDict
        );


    //- Destructor
    ~GAMGAgglomeration();


    // Member Functions

        // Access

            label size() const
            {
                return meshLevels_.size();
            }

            //- Return LDU mesh of given level
            const lduMesh& meshLevel(const label leveli) const;

            //- Return LDU interface addressing of given level
            const lduInterfacePtrsList& interfaceLevel
            (
                const label leveli
            ) const;

            //- Return cell restrict addressing of given level
            const labelField& restrictAddressing(const label leveli) const
            {
                return restrictAddressing_[leveli];
            }

            //- Return face restrict addressing of given level
            const labelList& faceRestrictAddressing(const label leveli) const
            {
                return faceRestrictAddressing_[leveli];
            }


        // Restriction and prolongation

            //- Restrict (integrate by summation) cell field
            template<class Type>
            void restrictField
            (
                Field<Type>& cf,
                const Field<Type>& ff,
                const label fineLevelIndex
            ) const;

            //- Restrict (integrate by summation) face field
            template<class Type>
            void restrictFaceField
            (
                Field<Type>& cf,
                const Field<Type>& ff,
                const label fineLevelIndex
            ) const;

            //- Prolong (interpolate by injection) cell field
            template<class Type>
            void prolongField
            (
                Field<Type>& ff,
                const Field<Type>& cf,
                const label coarseLevelIndex
            ) const;


        // Procesor agglomeration

            //- Collect and combine processor meshes into allMesh. Pass in
            //  communicator for combined mesh.
            void gatherMeshes
            (
                const labelList& procIDs,
                const label allMeshComm
            ) const;

            //- Combined mesh
            const lduPrimitiveMesh& allMesh() const;

            //- Mapping from processor to coarsestAllMesh cells
            const labelList& cellOffsets() const;

            //- Mapping from processor to coarsestAllMesh face
            const labelListList& faceMap() const;

            //- Mapping from processor to coarsestAllMesh boundary
            const labelListList& boundaryMap() const;

            //- Mapping from processor to coarsestAllMesh boundary face
            const labelListListList& boundaryFaceMap() const;

            //- Per remote processor the local equivalent. TBD: hack. Can be
            //  deleted once lduMatrices get received and merged into single
            //  one without constructing local ones
            const PtrList<lduMesh>& otherMeshes() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GAMGAgglomerationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
