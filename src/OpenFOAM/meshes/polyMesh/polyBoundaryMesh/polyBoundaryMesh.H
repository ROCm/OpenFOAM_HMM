/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::polyBoundaryMesh

Description
    A polyBoundaryMesh is a polyPatch list with additional search methods
    and registered IO.

SourceFiles
    polyBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef polyBoundaryMesh_H
#define polyBoundaryMesh_H

#include "polyPatchList.H"
#include "regIOobject.H"
#include "labelPair.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class wordRe;
class wordRes;

Ostream& operator<<(Ostream& os, const polyBoundaryMesh& pbm);


/*---------------------------------------------------------------------------*\
                      Class polyBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class polyBoundaryMesh
:
    public polyPatchList,
    public regIOobject
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Demand-driven: list of patch ids per face.
        mutable autoPtr<labelList> patchIDPtr_;

        //- Demand-driven: list of patch ids per group
        mutable autoPtr<HashTable<labelList>> groupIDsPtr_;

        //- Demand-driven: edges of neighbouring patches
        mutable autoPtr<List<labelPairList>> neighbourEdgesPtr_;


    // Private Member Functions

        //- Calculate geometry for the patches (transformation tensors etc.)
        void calcGeometry();

        //- Some patches have inGroup entries
        bool hasGroupIDs() const;

        //- Calculate group name to patch ids lookup
        void calcGroupIDs() const;

        //- No copy construct
        polyBoundaryMesh(const polyBoundaryMesh&) = delete;

        //- No copy assignment
        void operator=(const polyBoundaryMesh&) = delete;


public:

    //- Declare friendship with polyMesh
    friend class polyMesh;


    //- Runtime type information
    TypeName("polyBoundaryMesh");


    // Constructors

        //- Read constructor given IOobject and a polyMesh reference
        //  Note point pointers are unset, only used in copying meshes
        polyBoundaryMesh
        (
            const IOobject& io,
            const polyMesh& mesh
        );

        //- Construct given size
        polyBoundaryMesh
        (
            const IOobject& io,
            const polyMesh& mesh,
            const label size
        );

        //- Construct given polyPatchList
        polyBoundaryMesh
        (
            const IOobject& io,
            const polyMesh& mesh,
            const polyPatchList& ppl
        );


    //- Destructor
    ~polyBoundaryMesh() = default;


        //- Clear geometry at this level and at patches
        void clearGeom();

        //- Clear addressing at this level and at patches
        void clearAddressing();


    // Member Functions

        //- Return the mesh reference
        const polyMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Return a list of faceCells for each patch
        UPtrList<const labelUList> faceCells() const;

        //- Per patch the edges on the neighbouring patch.
        //  Is for every external edge the neighbouring patch and
        //  neighbouring (external) patch edge label. Note that edge indices
        //  are offset by nInternalEdges to keep it as much as possible
        //  consistent with coupled patch addressing (where coupling is by
        //  local patch face index). Only valid for singly connected
        //  polyBoundaryMesh and not parallel
        const List<labelPairList>& neighbourEdges() const;

        //- The number of patches before the first processor patch.
        label nNonProcessor() const;

        //- Return a list of patch names
        wordList names() const;

        //- Return a list of patch types
        wordList types() const;

        //- Return a list of physical types
        wordList physicalTypes() const;

        //- Return a list of patch start face indices
        labelList patchStarts() const;

        //- Return a list of patch sizes
        labelList patchSizes() const;

        //- Return a list of patch ranges
        List<labelRange> patchRanges() const;

        //- The start label of the boundary faces in the polyMesh face list
        //  Same as mesh.nInternalFaces()
        label start() const;

        //- The number of boundary faces in the underlying mesh
        //  Same as mesh.nBoundaryFaces()
        label nFaces() const;

        //- The face range for all boundary faces
        //  Spans [nInternalFaces, nFaces) of the underlying mesh
        labelRange range() const;

        //- Return the range used for boundary faces on patchi.
        //  Always returns an empty range for negative values of patchi,
        //  which allows safe use with findIndex or findPatchID.
        labelRange range(const label patchi) const;


        //- Return (sorted) patch indices for all matches.
        //  Optionally matches patch groups.
        //  A no-op (returns empty list) for an empty matcher
        labelList indices
        (
            const wordRe& matcher,
            const bool useGroups = true
        ) const;

        //- Return (sorted) patch indices for all matches.
        //  Optionally matches patch groups.
        //  A no-op (returns empty list) for an empty matcher
        labelList indices
        (
            const wordRes& matcher,
            const bool useGroups = true
        ) const;

        //- Return patch index for the first match, return -1 if not found
        //  A no-op (returns -1) for an empty key
        label findIndex(const wordRe& key) const;

        //- Find patch index given a name, return -1 if not found
        //  A no-op (returns -1) for an empty patchName
        label findPatchID
        (
            const word& patchName,
            const bool allowNotFound = true
        ) const;

        //- Find patch indices for a given polyPatch type
        template<class Type>
        labelHashSet findPatchIDs() const;

        //- Return patch index for a given face label
        label whichPatch(const label faceIndex) const;

        //- Per boundary face label the patch index
        const labelList& patchID() const;

        //- The patch indices per patch group
        const HashTable<labelList>& groupPatchIDs() const;

        //- Set/add group with patches
        void setGroup(const word& groupName, const labelUList& patchIDs);

        //- Return the set of patch IDs corresponding to the given names
        //  By default warns if given names are not found.
        //  Optionally matches to patchGroups as well as patchNames.
        labelHashSet patchSet
        (
            const UList<wordRe>& patchNames,
            const bool warnNotFound = true,
            const bool useGroups = true
        ) const;

        //- Match the patches to groups.
        //  Returns all the (fully matched) groups and any remaining
        //  unmatched patches.
        void matchGroups
        (
            const labelUList& patchIDs,
            wordList& groups,
            labelHashSet& nonGroupPatches
        ) const;

        //- Check whether all procs have all patches and in same order. Return
        //  true if in error.
        bool checkParallelSync(const bool report = false) const;

        //- Check boundary definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Correct polyBoundaryMesh after moving points
        void movePoints(const pointField& p);

        //- Correct polyBoundaryMesh after topology update
        void updateMesh();

        //- Reorders patches. Ordering does not have to be done in
        //  ascending or descending order. Reordering has to be unique.
        //  (is shuffle) If validBoundary calls updateMesh()
        //  after reordering to recalculate data (so call needs to be parallel
        //  sync in that case)
        void reorder(const labelUList& oldToNew, const bool validBoundary);

        //- writeData member function required by regIOobject
        virtual bool writeData(Ostream& os) const;

        //- Write using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool valid
        ) const;


    // Member Operators

        //- Return const and non-const reference to polyPatch by index.
        using polyPatchList::operator[];

        //- Return const reference to polyPatch by name.
        const polyPatch& operator[](const word& patchName) const;

        //- Return reference to polyPatch by name.
        polyPatch& operator[](const word& patchName);


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const polyBoundaryMesh& pbm);


    // Housekeeping

        //- Identical to the indices() method (AUG-2018)
        FOAM_DEPRECATED_FOR(2018-08, "indices() method")
        labelList findIndices(const wordRe& key, bool useGroups=true) const
        {
            return this->indices(key, useGroups);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "polyBoundaryMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
