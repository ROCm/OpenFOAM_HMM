/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::InteractionLists

Description

    Builds direct interaction list, specifying which local (real)
    cells are potentially in range of each other.

    Builds referred interaction list, specifying which cells are
    required to provideinteractions across coupled patched (cyclic or
    processor).  Generates referred cells, and refers particles to the
    correct processor, applying the appropriate transform.

    Simultaneous communication and computation is possible using:

    @verbatim
    PstreamBuffers pBufs(Pstream::nonBlocking);
    il_.sendReferredParticles(cellOccupancy_, pBufs);
    // Do other things
    il_.receiveReferredParticles(pBufs);
    @endverbatim

    Requiring data:

    @verbatim
    List<DynamicList<typename CloudType::parcelType*> > cellOccupancy_;
    @endverbatim

SourceFiles
    InteractionListsI.H
    InteractionLists.C
    InteractionListsIO.C

\*---------------------------------------------------------------------------*/

#ifndef InteractionLists_H
#define InteractionLists_H

#include "polyMesh.H"
#include "globalIndexAndTransform.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "mapDistribute.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class InteractionLists Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class InteractionLists
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Dummy cloud to give to particles
        Cloud<ParticleType> cloud_;

        //- mapDistribute to exchange particles
        autoPtr<mapDistribute> mapPtr_;

        //- Storage and encoding/decoding for all possible transforms
        //  of the geometry
        globalIndexAndTransform globalTransforms_;

        //- Maximum distance over which interactions will be detected
        scalar maxDistance_;

        //- Direct interaction list
        labelListList dil_;

        //- Wall faces on this processor that are in interaction range
        //  of each each cell, data subordinate to dil, i.e. if a cell
        //  B is on cell A's dil, then all of cell B's wall faces are
        //  considered posible interactions for cell A.
        labelListList directWallFaces_;

        //- Referred interaction list - which real cells are to be
        //  supplied with particles from the referred particle
        //  container with the same index.
        labelListList ril_;

        //- Inverse addressing for referred cells, specifies which
        //  referred cells (indices of entries in the ril_ and
        //  referredParticles_ lists) interact with the real cell
        //  indexed in this container.
        labelListList rilInverse_;

        //- Which cells are to be sent via the map, and an index
        //  specifying how they should be transformed.
        List<labelPair> cellIndexAndTransformToDistribute_;

        //- Referred particle container
        List<IDLList<ParticleType> > referredParticles_;


    // Private Member Functions

        //- Find the other processors which have interaction range
        //  extended bound boxes in range
        void findExtendedProcBbsInRange
        (
            const treeBoundBox& procBb,
            const List<treeBoundBox>& allExtendedProcBbs,
            const globalIndexAndTransform& globalTransforms,
            List<treeBoundBox>& extendedProcBbsInRange,
            List<label>& extendedProcBbsTransformIndex,
            List<label>& extendedProcBbsOrigProc
        );

        //- Build the mapDistribute
        void buildMap(const List<label>& toProc);

        //- Fill the referredParticles container with particles that
        //  will be referred
        void prepareParticlesToRefer
        (
            const List<DynamicList<ParticleType*> >& cellOccupancy
        );

        //- Prepare particle to be referred
        void prepareParticleToBeReferred
        (
            ParticleType* particle,
            labelPair giat
        );

        //- Write the referredParticles out to visualise
        void writeReferredParticleCloud();

        //- Disallow default bitwise copy construct
        InteractionLists(const InteractionLists&);

        //- Disallow default bitwise assignment
        void operator=(const InteractionLists&);


public:

    // Constructors

        //- Construct and create all information from the mesh
        InteractionLists
        (
            const polyMesh& mesh,
            scalar maxDistance
        );

    // Destructor

        ~InteractionLists();


    // Member Functions

        //- Prepare and send referred particles, nonBlocking communication
        void sendReferredParticles
        (
            const List<DynamicList<ParticleType*> >& cellOccupancy,
            PstreamBuffers& pBufs
        );

        //- Receive referred particles
        void receiveReferredParticles(PstreamBuffers& pBufs);

        // Access

            //- Return access to the mesh
            inline const polyMesh& mesh() const;

            //- Return access to the mapDistribute
            inline const mapDistribute& map() const;

            //- Return access to the globalTransforms
            inline const globalIndexAndTransform& globalTransforms() const;

            //- Return access to the direct interaction list
            inline const labelListList& dil() const;

            //- Return access to the direct wall faces
            inline const labelListList& directWallFaces() const;

            //- Return access to the referred interaction list
            inline const labelListList& ril() const;

            //- Return access to the inverse referred interaction list
            inline const labelListList& rilInverse() const;

            //- Return access to the cellIndexAndTransformToDistribute list
            inline const List<labelPair>&
            cellIndexAndTransformToDistribute() const;

            //- Return access to the referred particle container
            inline const List<IDLList<ParticleType> >&
            referredParticles() const;

            //- Return non-const access to the referred particle container
            inline List<IDLList<ParticleType> >& referredParticles();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InteractionListsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "InteractionLists.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
