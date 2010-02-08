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
    Foam::globalPoints

Description
    Calculates points shared by more than two processor patches or cyclic
    patches.

    Is used in globalMeshData. (this info is needed for point/edge
    communication where processor swaps are not enough to exchange data)

    Works purely topological and using local communication only.
    Needs:
      - domain to be one single domain (i.e. all faces can be reached through
        face-cell walk).
      - patch face ordering to be ok
      - f[0] ordering on patch faces to be ok.

    Works by constructing equivalence lists for all the points on processor
    patches. These list are in globalIndex numbering (so consecutively numbered
    per processor)
    E.g.
    @verbatim
          ((7 93)(4 731)(3 114))
    @endverbatim

    means point 93 on proc7 is connected to point 731 on proc4 and 114 on proc3.
    It then gets the lowest numbered processor (the 'master') to request a
    sharedPoint label from processor0 and it redistributes this label back to
    the other processors in the equivalence list.

    Algorithm:
        - get meshPoints of all my points on processor patches and initialize
          equivalence lists to this.
     loop
        - send to all neighbours in relative form:
            - patchFace
            - index in face
        - receive and convert into meshPoints. Add to to my equivalence lists.
        - mark meshPoints for which information changed.
        - send data for these meshPoints again
     endloop until nothing changes

    At this point one will have complete point-point connectivity for all
    points on processor patches. Now

        - (optional)remove point equivalences of size 2. These are
          just normal points shared between two neighbouring procPatches.
        - collect on each processor points for which it is the master
        - request number of sharedPointLabels from the Pstream::master.

    This information gets redistributed to all processors in a similar way
    as that in which the equivalence lists were collected:

        - initialize the indices of shared points I am the master for
     loop
        - send my known sharedPoints + meshPoints to all neighbours
        - receive from all neighbour. Find which meshPoint on my processor
          the sharedpoint is connected to
        - mark indices for which information has changed
     endloop until nothing changes.


    Note: the data held is either mesh point labels (construct from mesh only)
    or patch point labels (construct from mesh and patch).

SourceFiles
    globalPoints.C

\*---------------------------------------------------------------------------*/

#ifndef globalPoints_H
#define globalPoints_H

#include "DynamicList.H"
#include "Map.H"
#include "primitivePatch.H"
#include "edgeList.H"
#include "globalIndex.H"
#include "indirectPrimitivePatch.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyBoundaryMesh;
class cyclicPolyPatch;
class polyPatch;
class coupledPolyPatch;

/*---------------------------------------------------------------------------*\
                           Class globalPoints Declaration
\*---------------------------------------------------------------------------*/

class globalPoints
{
    // Static data members

        //- Offset to add to points (in globalIndices) originating from
        //  collocated coupled points.
        static const label fromCollocated;

        //- Distance to check whether points/faces are collocated.
        static const scalar mergeDist;


    // Private data

        //- Mesh reference
        const polyMesh& mesh_;

        //- Global numbering of points
        globalIndex globalIndices_;

        //- Sum of points on processor patches (unfiltered, point on 2 patches
        //  counts as 2)
        const label nPatchPoints_;

        //- All points on boundaries and their corresponding connected points
        //  on other processors.
        DynamicList<labelList> procPoints_;

        //- Map from mesh (or patch) point to index in procPoints
        Map<label> meshToProcPoint_;

        //- Shared points used by this processor (= global point number)
        labelList sharedPointAddr_;

        //- My mesh(or patch) points corresponding to the shared points
        labelList sharedPointLabels_;

        //- Total number of shared points.
        label nGlobalPoints_;


    // Private Member Functions

        //- Is identity transform?
        static bool noTransform(const tensor&, const scalar mergeDist);

        //- Return per face collocated status
        static PackedBoolList collocatedFaces
        (
            const coupledPolyPatch&,
            const scalar mergeDist
        );

        //- Return per point collocated status
        static PackedBoolList collocatedPoints
        (
            const coupledPolyPatch&,
            const scalar mergeDist
        );

        // Wrappers around global point numbering to add collocated bit

            //- Convert into globalIndices and add collocated bit
            label toGlobal(const label, const bool isCollocated) const;

            //- Is collocated bit set
            bool isCollocated(const label globalI) const;

            //- Remove collocated bit
            label removeCollocated(const label globalI) const;

            //- (remove collocated bit and) check if originates from local proc
            bool isLocal(const label globalI) const;

            //- (remove collocated bit and) get originating processor
            label whichProcID(const label globalI) const;

            //- (remove collocated bit and) convert to local number on processor
            label toLocal(const label procI, const label globalI) const;

            //- (remove collocated bit and) convert to local number on
            //  Pstream::myProcNo
            label toLocal(const label globalI) const;



        //- Count all points on processorPatches. Is all points for which
        //  information is collected.
        static label countPatchPoints(const polyBoundaryMesh&);

        //- Add information about patchPointI in relative indices to send
        //  buffers (patchFaces, indexInFace etc.)
        static void addToSend
        (
            const primitivePatch&,
            const label patchPointI,
            const labelList&,
            DynamicList<label>& patchFaces,
            DynamicList<label>& indexInFace,
            DynamicList<labelList>& allInfo
        );

        //- Merge info from neighbour into my data
        static bool mergeInfo
        (
            const labelList& nbrInfo,
            labelList& myInfo
        );

        //- From mesh point to 'local point'. Is the mesh point itself
        //  if meshToPatchPoint is empty.
        static label meshToLocalPoint
        (
            const Map<label>& meshToPatchPoint,
            const label meshPointI
        );

        //- Opposite of meshToLocalPoint.
        static label localToMeshPoint
        (
            const labelList& patchToMeshPoint,
            const label localPointI
        );

        //- Store (and merge) info for meshPointI
        bool storeInfo
        (
            const labelList& nbrInfo,
            const label localPointI,
            const bool isCollocated
        );

        void printProcPoints
        (
            const labelList& patchToMeshPoint,
            const labelList& pointInfo,
            Ostream& os
        ) const;

        //- Initialize procPoints_ to my patch points. allPoints = true:
        //  seed with all patch points, = false: only boundaryPoints().
        void initOwnPoints
        (
            const Map<label>& meshToPatchPoint,
            const bool allPoints,
            labelHashSet& changedPoints
        );

        //- Send subset of procPoints to neighbours
        void sendPatchPoints
        (
            const bool mergeSeparated,
            const Map<label>&,
            PstreamBuffers&,
            const labelHashSet&
        ) const;

        //- Receive neighbour points and merge into my procPoints.
        void receivePatchPoints
        (
            const bool mergeSeparated,
            const Map<label>&,
            PstreamBuffers&,
            labelHashSet&
        );

        //- Remove entries of size 2 where meshPoint is in provided Map.
        //  Used to remove normal face-face connected points.
        void remove(const labelList& patchToMeshPoint, const Map<label>&);

        //- Compact out unused elements of procPoints.
        void compact(const labelList& patchToMeshPoint);

        //- Get indices of point for which I am master (lowest numbered proc)
        labelList getMasterPoints(const labelList& patchToMeshPoint) const;


        //- Send subset of shared points to neighbours
        void sendSharedPoints
        (
            const bool mergeSeparated,
            PstreamBuffers&,
            const DynamicList<label>&
        ) const;

        //- Take over any local shared points
        void extendSharedPoints(const Map<label>&, DynamicList<label>&);

        //- Receive shared points and update subset.
        void receiveSharedPoints
        (
            const bool mergeSeparated,
            const Map<label>& meshToPatchPoint,
            const Map<label>& meshToShared,
            PstreamBuffers&,
            DynamicList<label>&
        );

        //- Should move into cyclicPolyPatch ordering problem
        //  keeps on giving problems.
        static edgeList coupledPoints(const cyclicPolyPatch&);

        //- Do all calculations.
        void calculateSharedPoints
        (
            const Map<label>&,
            const labelList&,
            const bool keepAllPoints,
            const bool mergeSeparated
        );

        //- Disallow default bitwise copy construct
        globalPoints(const globalPoints&);

        //- Disallow default bitwise assignment
        void operator=(const globalPoints&);


public:

        //- Declare name of the class and its debug switch
        ClassName("globalPoints");


    // Constructors

        //- Construct from mesh.
        //  keepAllPoints = false : filter out points that are on two
        //  neighbouring coupled patches (so can be swapped)
        //  mergeSeparated:
        //      true  : merge coupled points across separated patches.
        //      false : do not merge across coupled separated patches.
        globalPoints
        (
            const polyMesh& mesh,
            const bool keepAllPoints,
            const bool mergeSeparated
        );

        //- Construct from mesh and patch of coupled faces. Difference with
        //  construct from mesh only is that this stores the meshToProcPoint,
        //  procPoints and sharedPointLabels as patch local point labels
        //  instead of mesh point labels.
        globalPoints
        (
            const polyMesh& mesh,
            const indirectPrimitivePatch& coupledPatch,
            const bool keepAllPoints,
            const bool mergeSeparated
        );


    // Member Functions

        // Access

            //- From (mesh or patch) point to index in procPoints
            const Map<label>& meshToProcPoint() const
            {
                return meshToProcPoint_;
            }

            //- procPoints is per point the connected points (in global
            //  point numbers)
            const DynamicList<labelList>& procPoints() const
            {
                return procPoints_;
            }

            //- Global numbering of (mesh or patch) points
            const globalIndex& globalIndices() const
            {
                return globalIndices_;
            }

            //- shared points used by this processor (= global point number)
            const labelList& sharedPointAddr() const
            {
                return sharedPointAddr_;
            }

            //- my (mesh or patch)points corresponding to the shared points
            const labelList& sharedPointLabels() const
            {
                return sharedPointLabels_;
            }

            //- total number of shared points
            label nGlobalPoints() const
            {
                return nGlobalPoints_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
