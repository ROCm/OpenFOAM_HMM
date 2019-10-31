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
    Foam::ensightMesh

Description
    Encapsulation of volume meshes for writing in ensight format.

SourceFiles
    ensightMesh.C
    ensightMeshIO.C
    ensightMeshOptions.C

\*---------------------------------------------------------------------------*/

#ifndef ensightMesh_H
#define ensightMesh_H

#include "ensightCells.H"
#include "ensightFaces.H"
#include "ensightGeoFile.H"
#include "cellList.H"
#include "faceList.H"
#include "cellShapeList.H"
#include "HashTable.H"
#include "Map.H"
#include "scalarField.H"
#include "wordRes.H"
#include "globalIndex.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;
class ensightMesh;

/*---------------------------------------------------------------------------*\
                         Class ensightMesh Declaration
\*---------------------------------------------------------------------------*/

class ensightMesh
{
public:

    // Forward declarations
    class options;


private:

    // Private data

        //- Writer options
        const options* options_;

        //- Reference to the OpenFOAM mesh
        const fvMesh& mesh_;

        //- The volume cells (internalMesh)
        ensightCells meshCells_;

        //- Face elements per patch
        HashTable<ensightFaces> boundaryPatchFaces_;

        //- Face elements per faceZone
        HashTable<ensightFaces> faceZoneFaces_;

        //- The list of patches to be output
        Map<word> patchLookup_;

        //- Track if it needs an update
        mutable bool needsUpdate_;


        // Parallel merged points

        //- Global numbering for merged points
        autoPtr<globalIndex> globalPointsPtr_;

        //- From mesh point to global merged point
        labelList pointToGlobal_;

        //- Local points that are unique
        labelList uniquePointMap_;


    // Private Member Functions

        //- Clear some storage
        void clear();


        //- Inplace renumber of cell-shapes
        static cellShapeList& renumberShapes
        (
            cellShapeList& shapes,
            const labelUList& pointToGlobal
        );

        //- Copy and return renumbered cell-shapes
        static cellShapeList renumberShapes
        (
            const cellShapeList& shapes,
            const labelUList& addr,
            const labelUList& pointToGlobal
        );

        //- Write list of faces
        static void writeFaceList
        (
            const faceList& faces,
            ensightGeoFile& os
        );

        //- Write list of faces
        static void writeFaceList
        (
            const UIndirectList<face>& faces,
            ensightGeoFile& os
        );

        //- Return sizes of faces in the list
        static labelList getFaceSizes
        (
            const faceList& faces
        );

        //- Return sizes of faces in the list
        static labelList getFaceSizes
        (
            const UIndirectList<face>& faces
        );

        //- Write sizes of faces in the list
        static void writeFaceSizes
        (
            const faceList& faces,
            ensightGeoFile& os
        );

        //- Write sizes of faces in the list
        static void writeFaceSizes
        (
            const UIndirectList<face>& faces,
            ensightGeoFile& os
        );

        //- Write cell connectivity via cell shapes
        static void writeCellShapes
        (
            const cellShapeList& shapes,
            ensightGeoFile& os
        );

        //- Return the number of faces per poly element
        static labelList getPolysNFaces
        (
            const labelUList& polys,
            const cellList& cellFaces
        );

        //- Write the number of faces per poly element
        static void writePolysNFaces
        (
            const labelUList& polys,
            const cellList& cellFaces,
            ensightGeoFile& os
        );

        //- Return the number of points per poly element
        static labelList getPolysNPointsPerFace
        (
            const labelUList& polys,
            const cellList& cellFaces,
            const faceList& faces
        );

        //- Write the number of points per poly element
        static void writePolysNPointsPerFace
        (
            const labelUList& polys,
            const cellList& cellFaces,
            const faceList& faces,
            ensightGeoFile& os
        );

        //- Write the point ids per poly element
        static void writePolysPoints
        (
            const labelUList& addr,
            const cellList& cellFaces,
            const faceList& faces,
            const labelList& faceOwner,
            ensightGeoFile& os
        );

        //- Write the poly connectivity
        void writePolysConnectivity
        (
            const labelUList& polys,
            const labelList& pointToGlobal,
            ensightGeoFile&
        ) const;

        //- Write the regular cell connectivity for all types
        void writeCellConnectivity
        (
            const ensightCells& ensCells,
            const labelList& pointToGlobal,
            ensightGeoFile& os
        ) const;

        //- Write the regular cell connectivity for specified type
        void writeCellConnectivity
        (
            ensightCells::elemType elemType,
            const ensightCells& ensCells,
            const labelList& pointToGlobal,
            ensightGeoFile& os
        ) const;

        //- Write the regular face connectivity for specified type and
        //- and specified faces
        void writeFaceConnectivity
        (
            ensightFaces::elemType elemType,
            const label nTotal,
            const faceList& faces,
            const labelUList& addr,
            ensightGeoFile&
        ) const;

        //- Write the regular face connectivity for specified type
        void writeFaceConnectivity
        (
            ensightFaces::elemType elemType,
            const label nTotal,
            const faceList& faces,
            ensightGeoFile& os
        ) const;


        void writeFaceConnectivity
        (
            const ensightFaces& ensFaces,
            const faceList& faces,
            ensightGeoFile& os,
            const bool raw = false
        ) const;


        void writeAllPoints
        (
            const label partId,
            const word& ensightPartName,
            const label nTotal,
            const pointField& uniquePoints,
            ensightGeoFile&
        ) const;


        //- No copy construct
        ensightMesh(const ensightMesh&) = delete;

        //- No copy assignment
        void operator=(const ensightMesh&) = delete;


public:

    // Constructors

        //- Construct from components
        ensightMesh(const fvMesh& mesh, const options& opts);

        //- Construct from fvMesh with all default options, binary output
        explicit ensightMesh(const fvMesh& mesh);

        //- Construct from fvMesh with all default options and specified format
        ensightMesh(const fvMesh& mesh, const IOstream::streamFormat format);


    //- Destructor
    ~ensightMesh();


    // Member Functions

    // Access

        //- Reference to the underlying fvMesh
        inline const fvMesh& mesh() const;

        //- Reference to the writer/mesh options
        inline const ensightMesh::options& option() const;

        //- Ascii/Binary file output
        inline IOstream::streamFormat format() const;

        //- Using internal?
        inline bool useInternalMesh() const;

        //- Using boundary?
        inline bool useBoundaryMesh() const;

        //- The volume cells (internalMesh)
        inline const ensightCells& meshCells() const;

        //- The list of patches to be output
        inline const Map<word>& patches() const;

        //- Face elements per selected patch
        inline const HashTable<ensightFaces>& boundaryPatchFaces() const;

        //- Face elements per selected faceZone.
        //  To be output in sorted order.
        inline const HashTable<ensightFaces>& faceZoneFaces() const;


    // Parallel point merging

        //- Global numbering for merged points
        const globalIndex& globalPoints() const
        {
            return globalPointsPtr_();
        }

        //- From mesh point to global merged point
        const labelList& pointToGlobal() const
        {
            return pointToGlobal_;
        }

        //- Local points that are unique
        const labelList& uniquePointMap() const
        {
            return uniquePointMap_;
        }


    // Other

        //- Does the content need an update?
        bool needsUpdate() const;

        //- Mark as needing an update.
        //  May also free up unneeded data.
        //  Return false if already marked as expired.
        bool expire();

        //- Update for new mesh
        void correct();


    // Output

        //- Write to file
        inline void write(autoPtr<ensightGeoFile>& os) const;

        //- Write to file
        void write(ensightGeoFile& os) const;

};


//- Configuration options for the ensightMesh
class ensightMesh::options
{
    //- Ascii/Binary file output
    IOstream::streamFormat format_;

    //- Create in 'expired' mode
    bool lazy_;

    //- Use the internal mesh
    bool internal_;

    //- Use the boundary mesh
    bool boundary_;

    //- Output of selected patches only
    wordRes patchPatterns_;

    //- Output of selected faceZones
    wordRes faceZonePatterns_;


public:

    // Constructors

        //- Construct for binary output
        options();

        //- Construct for specified format
        explicit options(IOstream::streamFormat format);


    // Member Functions

    // Access

        //- File output format (ascii | binary)
        IOstream::streamFormat format() const;

        //- Lazy creation? (ie, ensightMesh starts as needsUpdate)
        bool lazy() const;

        //- Using internal?
        bool useInternalMesh() const;

        //- Using boundary?
        bool useBoundaryMesh() const;

        //- Using faceZones?
        bool useFaceZones() const;

        //- Selection of patches. Empty if unspecified.
        const wordRes& patchSelection() const;

        //- Selection of faceZones. Empty if unspecified.
        const wordRes& faceZoneSelection() const;


    // Edit

        //- Reset to defaults
        void reset();

        //- Lazy creation - ensightMesh starts as needsUpdate.
        void lazy(bool beLazy);

        //- Alter the useBoundaryMesh state
        void useInternalMesh(bool on);

        //- Alter the useBoundaryMesh state
        void useBoundaryMesh(bool on);

        //- Define patch selection matcher
        void patchSelection(const UList<wordRe>& patterns);

        //- Define patch selection matcher
        void patchSelection(List<wordRe>&& patterns);

        //- Define faceZone selection matcher
        void faceZoneSelection(const UList<wordRe>& patterns);

        //- Define faceZone selection matcher
        void faceZoneSelection(List<wordRe>&& patterns);


    // Housekeeping

        //- Older name for useBoundaryMesh()
        //  \deprecated OCT-2018
        bool usePatches() const { return useBoundaryMesh(); }

        //- Older name for useBoundaryMesh()
        //  \deprecated OCT-2018
        void noPatches(bool off) { useBoundaryMesh(!off); }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightMeshI.H"

#endif

// ************************************************************************* //
