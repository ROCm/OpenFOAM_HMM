/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Foam::ensightMesh

Description

SourceFiles
    ensightMesh.C

\*---------------------------------------------------------------------------*/

#ifndef ensightMesh_H
#define ensightMesh_H

#include "cellSets.H"
#include "faceSets.H"
#include "HashTable.H"
#include "HashSet.H"
#include "fvMesh.H"
#include "OFstream.H"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class argList;

/*---------------------------------------------------------------------------*\
                           Class ensightMesh Declaration
\*---------------------------------------------------------------------------*/

class ensightMesh
{
    // Private data

        //- Set binary file output
        bool binary_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ensightMesh(const ensightMesh&);

        //- Disallow default bitwise assignment
        void operator=(const ensightMesh&);

        void writePoints
        (
            const scalarField& pointsComponent,
            OFstream& ensightGeometryFile
        ) const;

        cellShapeList map
        (
            const cellShapeList& cellShapes,
            const labelList& prims
        ) const;

        cellShapeList map
        (
            const cellShapeList& cellShapes,
            const labelList& hexes,
            const labelList& wedges
        ) const;

        void writePrims
        (
            const cellShapeList& cellShapes,
            const label pointOffset,
            OFstream& ensightGeometryFile
        ) const;

        void writePolys
        (
            const labelList& polys,
            const cellList& cellFaces,
            const faceList& faces,
            const label pointOffset,
            OFstream& ensightGeometryFile
        ) const;

        void writeAllPrims
        (
            const char* key,
            const label nPrims,
            const cellShapeList& cellShapes,
            const labelList& pointOffsets,
            OFstream& ensightGeometryFile
        ) const;

        void writeFacePrims
        (
            const char* key,
            const faceList& patchFaces,
            const label pointOffset,
            OFstream& ensightGeometryFile
        ) const;

        faceList map
        (
            const faceList& patchFaces,
            const labelList& prims
        ) const;

        void writeAllFacePrims
        (
            const char* key,
            const labelList& prims,
            const label nPrims,
            const faceList& patchFaces,
            const labelList& pointOffsets,
            const labelList& patchProcessors,
            OFstream& ensightGeometryFile
        ) const;

        void writeAscii
        (
            const fileName& postProcPath,
            const word& prepend,
            const label timeIndex,
            Ostream& ensightCaseFile
        ) const;

        void writeBinary
        (
            const fileName& postProcPath,
            const word& prepend,
            const label timeIndex,
            Ostream& ensightCaseFile
        ) const;

        void writePrimsBinary
        (
            const cellShapeList& cellShapes,
            const label pointOffset,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeAllPrimsBinary
        (
            const char* key,
            const label nPrims,
            const cellShapeList& cellShapes,
            const labelList& pointOffsets,
            std::ofstream& ensightGeometryFile
        ) const;

        void writePolysBinary
        (
            const labelList& polys,
            const cellList& cellFaces,
            const faceList& faces,
            const label pointOffset,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeAllFacePrimsBinary
        (
            const char* key,
            const labelList& prims,
            const label nPrims,
            const faceList& patchFaces,
            const labelList& pointOffsets,
            const labelList& patchProcessors,
            std::ofstream& ensightGeometryFile
        ) const;

        void writeFacePrimsBinary
        (
            const char* key,
            const faceList& patchFaces,
            const label pointOffset,
            std::ofstream& ensightGeometryFile
        ) const;


public:

    // Public data

        const fvMesh& mesh;
        cellSets meshCellSets;
        List<faceSets> boundaryFaceSets;

        HashTable<labelList> allPatchNames;
        HashTable<label> patchIndices;
        wordHashSet patchNames;

        class nFacePrims
        {
        public:

            label nPoints;
            label nTris;
            label nQuads;
            label nPolys;

            nFacePrims()
            :
                nPoints(0),
                nTris(0),
                nQuads(0),
                nPolys(0)
            {}
        };

        HashTable<nFacePrims> nPatchPrims;


    // Constructors

        //- Construct from fvMesh
        ensightMesh(const fvMesh&, const argList& args, const bool binary);


    // Destructor

        ~ensightMesh();


    // Member Functions

        void write
        (
            const fileName& postProcPath,
            const word& prepend,
            const label timeIndex,
            Ostream& ensightCaseFile
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
