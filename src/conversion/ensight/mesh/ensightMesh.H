/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
#include "wordReList.H"
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
            cellShapeList&,
            const labelUList& pointToGlobal
        );

        static cellShapeList map
        (
            const cellShapeList&,
            const labelUList& prims,
            const labelUList& pointToGlobal
        );

        //- Write list of faces
        static void writeFaceList
        (
            const faceList&,
            ensightGeoFile&
        );

        //- Write list of faces
        static void writeFaceList
        (
            const UIndirectList<face>&,
            ensightGeoFile&
        );

        //- Write sizes of faces in the list
        static void writeFaceSizes
        (
            const faceList&,
            ensightGeoFile&
        );

        //- Write sizes of faces in the list
        static void writeFaceSizes
        (
            const UIndirectList<face>&,
            ensightGeoFile&
        );

        //- Write cell connectivity via shell shapes
        static void writeCellShapes
        (
            const cellShapeList&,
            ensightGeoFile&
        );

        void writePolysNFaces
        (
            const labelList& polys,
            const cellList& cellFaces,
            ensightGeoFile&
        ) const;

        void writePolysNPointsPerFace
        (
            const labelList& polys,
            const cellList& cellFaces,
            const faceList& faces,
            ensightGeoFile&
        ) const;

        void writePolysPoints
        (
            const labelList& polys,
            const cellList& cellFaces,
            const faceList& faces,
            const labelList& faceOwner,
            ensightGeoFile&
        ) const;

        void writePolysConnectivity
        (
            const labelList& addr,
            const labelList& pointToGlobal,
            ensightGeoFile&
        ) const;

        void writeCellConnectivity
        (
            const ensightCells&,
            const labelList& pointToGlobal,
            ensightGeoFile&
        ) const;

        void writeCellConnectivity
        (
            ensightCells::elemType elemType,
            const ensightCells&,
            const labelList& pointToGlobal,
            ensightGeoFile&
        ) const;

        void writeFaceConnectivity
        (
            ensightFaces::elemType elemType,
            const label nTotal,
            const faceList& faceLst,
            const labelList& addr,
            ensightGeoFile&
        ) const;


        void writeFaceConnectivity
        (
            ensightFaces::elemType elemType,
            const label nTotal,
            const faceList& faceLst,
            ensightGeoFile&
        ) const;


        void writeFaceConnectivity
        (
            const ensightFaces&,
            const faceList& faceLst,
            ensightGeoFile&,
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


        //- Disallow default bitwise copy construct
        ensightMesh(const ensightMesh&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const ensightMesh&) = delete;


public:

    // Constructors

        //- Construct from components
        ensightMesh
        (
            const fvMesh& mesh,
            const options& opts
        );

        //- Construct from fvMesh with all default options
        ensightMesh
        (
            const fvMesh& mesh,
            const IOstream::streamFormat format = IOstream::BINARY
        );


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

        //- Using internalMesh?
        inline bool useInternalMesh() const;

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


    // I-O

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

    //- Suppress patches
    bool noPatches_;

    //- Output selected patches only
    autoPtr<wordReList> patchPatterns_;

    //- Output selected faceZones
    autoPtr<wordReList> faceZonePatterns_;

public:

    // Constructors

        //- Construct with the specified format (default is binary)
        options(IOstream::streamFormat format = IOstream::BINARY);


    // Member Functions

    // Access

        //- Ascii/Binary file output
        IOstream::streamFormat format() const;

        //- Lazy creation? (ie, ensightMesh starts as needsUpdate)
        bool lazy() const;

        //- Using internalMesh?
        bool useInternalMesh() const;

        //- Using patches?
        bool usePatches() const;

        //- Using faceZones?
        bool useFaceZones() const;

        //- Selection of patches in effect?
        bool usePatchSelection() const;

        //- Selection of patches - null reference if not available
        const wordReList& patchSelection() const;

        //- Selection of faceZones - null reference if not available
        const wordReList& faceZoneSelection() const;


    // Edit

        //- Reset to defaults
        void reset();

        //- Lazy creation - ensightMesh starts as needsUpdate.
        void lazy(const bool);

        //- Alter the patches/no-patches state
        void noPatches(const bool);

        //- Define patch selection matcher
        void patchSelection(const wordReList&);

        //- Define faceZone selection matcher
        void faceZoneSelection(const wordReList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightMeshI.H"

#endif

// ************************************************************************* //
