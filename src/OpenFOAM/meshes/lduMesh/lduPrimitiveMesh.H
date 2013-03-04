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
    Foam::lduPrimitiveMesh

Description
    Simplest contrete lduMesh which stores the addressing needed by lduMatrix.

SourceFiles
    lduPrimitiveMesh.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveMesh_H
#define lduPrimitiveMesh_H

#include "lduMesh.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class lduPrimitiveMesh Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveMesh
:
    public lduMesh,
    public lduAddressing
{
    // Private data

        //- Lower addressing
        labelList lowerAddr_;

        //- Upper addressing
        labelList upperAddr_;

        //- Patch to internal addressing
        labelListList patchAddr_;

        //- List of pointers for each patch
        //  with only those pointing to interfaces being set
        lduInterfacePtrsList interfaces_;

        //- Concrete interfaces if above interfaces were created locally
        PtrList<lduInterface> primitiveInterfaces_;

        //- Patch field evaluation schedule
        lduSchedule patchSchedule_;

        //- Communicator to use for any parallel communication
        const label comm_;

    // Private Member Functions

        //- Get size of all meshes
        static label size(const PtrList<lduMesh>&);

        //- Check if in upper-triangular ordering
        static void checkUpperTriangular
        (
            const label size,
            const labelUList& l,
            const labelUList& u
        );

        //- Get non-scheduled send/receive schedule
        template<class ProcPatch>
        static lduSchedule nonBlockingSchedule(const lduInterfacePtrsList&);

        //- Disallow default bitwise copy construct
        lduPrimitiveMesh(const lduPrimitiveMesh&);

        //- Disallow default bitwise assignment
        void operator=(const lduPrimitiveMesh&);


public:

    // Constructors

        //- Construct from components as copies
        lduPrimitiveMesh
        (
            const label nCells,
            const labelUList& l,
            const labelUList& u,
            const labelListList& pa,
            const lduInterfacePtrsList interfaces,
            const lduSchedule& ps,
            const label comm
        );

        //- Construct from components and re-use storage as specified.
        lduPrimitiveMesh
        (
            const label nCells,
            labelList& l,
            labelList& u,
            labelListList& pa,
            lduInterfacePtrsList interfaces,
            const lduSchedule& ps,
            const label comm,
            bool reUse
        );

        //- Construct by combining multiple meshes. The meshes come from
        //  processors procIDs:
        //  procIDs[0] : local processor (myMesh)
        //  procIDs[i] : processor where otherMeshes[i-1] comes from
        //  procAgglomMap : for every processor which processor it agglomerates
        //                  onto. The new processor numbers are in compact
        //                  numbering (so ranks in communicator comm), i.e.
        //                  similar to cell-restrict-addressing.
        //                  We need this information to be able to map
        //                  inter-processor interfaces
        lduPrimitiveMesh
        (
            const label comm,
            const labelList& procAgglomMap,

            const labelList& procIDs,               // procIDs[0] = myMesh
            const lduMesh& myMesh,
            const PtrList<lduMesh>& otherMeshes,

            labelList& cellOffsets,
            labelListList& faceMap,
            labelListList& boundaryMap,
            labelListListList& boundaryFaceMap
        );


    //- Destructor
    virtual ~lduPrimitiveMesh()
    {}


    // Member Functions

        // Access

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const
            {
                return *this;
            }

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return interfaces_;
            }

            //- Return communicator used for parallel communication
            virtual label comm() const
            {
                return comm_;
            }

            //- Return Lower addressing
            virtual const labelUList& lowerAddr() const
            {
                return lowerAddr_;
            }

            //- Return Upper addressing
            virtual const labelUList& upperAddr() const
            {
                return upperAddr_;
            }

            //- Return patch addressing
            virtual const labelUList& patchAddr(const label i) const
            {
                return patchAddr_[i];
            }

            //- Return patch evaluation schedule
            virtual const lduSchedule& patchSchedule() const
            {
                return patchSchedule_;
            }


        // Helper

            //- Select either mesh0 (meshI is 0) or otherMeshes[meshI-1]
            static const lduMesh& mesh
            (
                const lduMesh& mesh0,
                const PtrList<lduMesh>& otherMeshes,
                const label meshI
            );

            //- Gather meshes from other processors onto procIDs[0].
            //  Hardcoded for GAMGInterface.
            static void gather
            (
                const lduMesh& mesh,
                const labelList& procIDs,
                PtrList<lduMesh>& otherMeshes
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "lduPrimitiveMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
