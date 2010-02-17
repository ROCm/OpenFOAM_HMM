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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "ensightMesh.H"
#include "fvMesh.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "cellModeller.H"
#include "IOmanip.H"
#include "itoa.H"
#include "ensightWriteBinary.H"
#include "mapDistribute.H"

#include <fstream>

// * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightMesh::ensightMesh
(
    const fvMesh& mesh,
    const argList& args,
    const bool binary
)
:
    mesh_(mesh),
    binary_(binary),
    patchPartOffset_(2),
    meshCellSets_(mesh_.nCells()),
    boundaryFaceSets_(mesh_.boundary().size()),
    allPatchNames_(0),
    patchNames_(0),
    nPatchPrims_(0)
{
    const cellShapeList& cellShapes = mesh.cellShapes();

    const cellModel& tet = *(cellModeller::lookup("tet"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& wedge = *(cellModeller::lookup("wedge"));
    const cellModel& hex = *(cellModeller::lookup("hex"));

    if (!args.optionFound("noPatches"))
    {
        // Patches are output. Check that they're synced.
        mesh_.boundaryMesh().checkParallelSync(true);

        allPatchNames_ = wordList::subList
        (
            mesh_.boundaryMesh().names(),
            mesh_.boundary().size()
          - mesh_.globalData().processorPatches().size()
        );

        if (args.optionFound("patches"))
        {
            wordList patchNameList(args.optionLookup("patches")());

            if (patchNameList.empty())
            {
                patchNameList = allPatchNames_;
            }

            forAll(patchNameList, i)
            {
                patchNames_.insert(patchNameList[i]);
            }
        }
    }

    if (patchNames_.size())
    {
        // no internalMesh
        patchPartOffset_ = 1;
    }
    else
    {
        // Count the shapes
        labelList& tets = meshCellSets_.tets;
        labelList& pyrs = meshCellSets_.pyrs;
        labelList& prisms = meshCellSets_.prisms;
        labelList& wedges = meshCellSets_.wedges;
        labelList& hexes = meshCellSets_.hexes;
        labelList& polys = meshCellSets_.polys;

        label nTets = 0;
        label nPyrs = 0;
        label nPrisms = 0;
        label nWedges = 0;
        label nHexes = 0;
        label nPolys = 0;

        forAll(cellShapes, cellI)
        {
            const cellShape& cellShape = cellShapes[cellI];
            const cellModel& cellModel = cellShape.model();

            if (cellModel == tet)
            {
                tets[nTets++] = cellI;
            }
            else if (cellModel == pyr)
            {
                pyrs[nPyrs++] = cellI;
            }
            else if (cellModel == prism)
            {
                prisms[nPrisms++] = cellI;
            }
            else if (cellModel == wedge)
            {
                wedges[nWedges++] = cellI;
            }
            else if (cellModel == hex)
            {
                hexes[nHexes++] = cellI;
            }
            else
            {
                polys[nPolys++] = cellI;
            }
        }

        tets.setSize(nTets);
        pyrs.setSize(nPyrs);
        prisms.setSize(nPrisms);
        wedges.setSize(nWedges);
        hexes.setSize(nHexes);
        polys.setSize(nPolys);

        meshCellSets_.nTets = nTets;
        reduce(meshCellSets_.nTets, sumOp<label>());

        meshCellSets_.nPyrs = nPyrs;
        reduce(meshCellSets_.nPyrs, sumOp<label>());

        meshCellSets_.nPrisms = nPrisms;
        reduce(meshCellSets_.nPrisms, sumOp<label>());

        meshCellSets_.nHexesWedges = nHexes + nWedges;
        reduce(meshCellSets_.nHexesWedges, sumOp<label>());

        meshCellSets_.nPolys = nPolys;
        reduce(meshCellSets_.nPolys, sumOp<label>());
    }

    if (!args.optionFound("noPatches"))
    {
        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].size())
            {
                const polyPatch& p = mesh.boundaryMesh()[patchi];

                labelList& tris = boundaryFaceSets_[patchi].tris;
                labelList& quads = boundaryFaceSets_[patchi].quads;
                labelList& polys = boundaryFaceSets_[patchi].polys;

                tris.setSize(p.size());
                quads.setSize(p.size());
                polys.setSize(p.size());

                label nTris = 0;
                label nQuads = 0;
                label nPolys = 0;

                forAll(p, faceI)
                {
                    const face& f = p[faceI];

                    if (f.size() == 3)
                    {
                        tris[nTris++] = faceI;
                    }
                    else if (f.size() == 4)
                    {
                        quads[nQuads++] = faceI;
                    }
                    else
                    {
                        polys[nPolys++] = faceI;
                    }
                }

                tris.setSize(nTris);
                quads.setSize(nQuads);
                polys.setSize(nPolys);
            }
        }
    }


    forAll(allPatchNames_, patchi)
    {
        const word& patchName = allPatchNames_[patchi];
        nFacePrimitives nfp;

        if (patchNames_.empty() || patchNames_.found(patchName))
        {
            if (mesh.boundary()[patchi].size())
            {
                nfp.nTris   = boundaryFaceSets_[patchi].tris.size();
                nfp.nQuads  = boundaryFaceSets_[patchi].quads.size();
                nfp.nPolys  = boundaryFaceSets_[patchi].polys.size();
            }
        }

        reduce(nfp.nTris, sumOp<label>());
        reduce(nfp.nQuads, sumOp<label>());
        reduce(nfp.nPolys, sumOp<label>());

        nPatchPrims_.insert(patchName, nfp);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightMesh::~ensightMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightMesh::writePoints
(
    const scalarField& pointsComponent,
    OFstream& ensightGeometryFile
) const
{
    forAll(pointsComponent, pointI)
    {
        ensightGeometryFile<< setw(12) << float(pointsComponent[pointI]) << nl;
    }
}


Foam::cellShapeList Foam::ensightMesh::map
(
    const cellShapeList& cellShapes,
    const labelList& prims,
    const labelList& pointToGlobal
) const
{
    cellShapeList mcsl(prims.size());

    forAll(prims, i)
    {
        mcsl[i] = cellShapes[prims[i]];
        inplaceRenumber(pointToGlobal, mcsl[i]);
    }

    return mcsl;
}


Foam::cellShapeList Foam::ensightMesh::map
(
    const cellShapeList& cellShapes,
    const labelList& hexes,
    const labelList& wedges,
    const labelList& pointToGlobal
) const
{
    cellShapeList mcsl(hexes.size() + wedges.size());

    forAll(hexes, i)
    {
        mcsl[i] = cellShapes[hexes[i]];
        inplaceRenumber(pointToGlobal, mcsl[i]);
    }

    label offset = hexes.size();

    const cellModel& hex = *(cellModeller::lookup("hex"));
    labelList hexLabels(8);

    forAll(wedges, i)
    {
        const cellShape& cellPoints = cellShapes[wedges[i]];

        hexLabels[0] = cellPoints[0];
        hexLabels[1] = cellPoints[1];
        hexLabels[2] = cellPoints[0];
        hexLabels[3] = cellPoints[2];
        hexLabels[4] = cellPoints[3];
        hexLabels[5] = cellPoints[4];
        hexLabels[6] = cellPoints[6];
        hexLabels[7] = cellPoints[5];

        mcsl[i + offset] = cellShape(hex, hexLabels);
        inplaceRenumber(pointToGlobal, mcsl[i + offset]);
    }

    return mcsl;
}


void Foam::ensightMesh::writePrims
(
    const cellShapeList& cellShapes,
    OFstream& ensightGeometryFile
) const
{
    forAll(cellShapes, i)
    {
        const cellShape& cellPoints = cellShapes[i];

        forAll(cellPoints, pointI)
        {
            ensightGeometryFile
                << setw(10)
                << cellPoints[pointI] + 1;
        }
        ensightGeometryFile << nl;
    }
}


void Foam::ensightMesh::writePrimsBinary
(
    const cellShapeList& cellShapes,
    std::ofstream& ensightGeometryFile
) const
{
    // Create a temp int array
    int numElem;

    numElem = cellShapes.size();

    if (cellShapes.size())
    {
        // All the cellShapes have the same number of elements!
        int numIntElem = cellShapes.size()*cellShapes[0].size();
        List<int> temp(numIntElem);

        int n = 0;

        forAll(cellShapes, i)
        {
            const cellShape& cellPoints = cellShapes[i];

            forAll(cellPoints, pointI)
            {
                temp[n] = cellPoints[pointI] + 1;
                n++;
            }
        }

        ensightGeometryFile.write
        (
            reinterpret_cast<char*>(temp.begin()),
            numIntElem*sizeof(int)
        );
    }
}


void Foam::ensightMesh::writePolysNFaces
(
    const labelList& polys,
    const cellList& cellFaces,
    OFstream& ensightGeometryFile
) const
{
    forAll(polys, i)
    {
        ensightGeometryFile
            << setw(10) << cellFaces[polys[i]].size() << nl;
    }
}


void Foam::ensightMesh::writePolysNPointsPerFace
(
    const labelList& polys,
    const cellList& cellFaces,
    const faceList& faces,
    OFstream& ensightGeometryFile
) const
{
    forAll(polys, i)
    {
        const labelList& cf = cellFaces[polys[i]];

        forAll(cf, faceI)
        {
            ensightGeometryFile
                << setw(10) << faces[cf[faceI]].size() << nl;
        }
    }
}


void Foam::ensightMesh::writePolysPoints
(
    const labelList& polys,
    const cellList& cellFaces,
    const faceList& faces,
    OFstream& ensightGeometryFile
) const
{
    forAll(polys, i)
    {
        const labelList& cf = cellFaces[polys[i]];

        forAll(cf, faceI)
        {
            const face& f = faces[cf[faceI]];

            forAll(f, pointI)
            {
                ensightGeometryFile << setw(10) << f[pointI] + 1;
            }
            ensightGeometryFile << nl;
        }
    }
}


void Foam::ensightMesh::writeAllPolys
(
    const labelList& pointToGlobal,
    OFstream& ensightGeometryFile
) const
{
    if (meshCellSets_.nPolys)
    {
        const cellList& cellFaces = mesh_.cells();
        // Renumber faces to use global point numbers
        faceList faces(mesh_.faces());
        forAll(faces, i)
        {
            inplaceRenumber(pointToGlobal, faces[i]);
        }

        if (Pstream::master())
        {
            ensightGeometryFile
                << "nfaced" << nl << setw(10) << meshCellSets_.nPolys << nl;
        }

        // Number of faces for each poly cell
        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                toMaster<< meshCellSets_.polys << cellFaces;
            }

            pBufs.finishedSends();

            if (Pstream::master())
            {
                // Master
                writePolysNFaces
                (
                    meshCellSets_.polys,
                    cellFaces,
                    ensightGeometryFile
                );
                // Slaves
                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    UIPstream fromSlave(slave, pBufs);
                    labelList polys(fromSlave);
                    cellList cellFaces(fromSlave);

                    writePolysNFaces
                    (
                        polys,
                        cellFaces,
                        ensightGeometryFile
                    );
                }
            }
        }


        // Number of points for each face of the above list
        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                toMaster<< meshCellSets_.polys << cellFaces << faces;
            }
            pBufs.finishedSends();

            if (Pstream::master())
            {
                // Master
                writePolysNPointsPerFace
                (
                    meshCellSets_.polys,
                    cellFaces,
                    faces,
                    ensightGeometryFile
                );
                // Slaves
                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    UIPstream fromSlave(slave, pBufs);
                    labelList polys(fromSlave);
                    cellList cellFaces(fromSlave);
                    faceList faces(fromSlave);

                    writePolysNPointsPerFace
                    (
                        polys,
                        cellFaces,
                        faces,
                        ensightGeometryFile
                    );
                }
            }
        }


        // List of points id for each face of the above list
        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                toMaster<< meshCellSets_.polys << cellFaces << faces;
            }

            pBufs.finishedSends();

            if (Pstream::master())
            {
                // Master
                writePolysPoints
                (
                    meshCellSets_.polys,
                    cellFaces,
                    faces,
                    ensightGeometryFile
                );
                // Slaves
                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    UIPstream fromSlave(slave, pBufs);
                    labelList polys(fromSlave);
                    cellList cellFaces(fromSlave);
                    faceList faces(fromSlave);

                    writePolysPoints
                    (
                        polys,
                        cellFaces,
                        faces,
                        ensightGeometryFile
                    );
                }
            }
        }
    }
}


void Foam::ensightMesh::writePolysNFacesBinary
(
    const labelList& polys,
    const cellList& cellFaces,
    std::ofstream& ensightGeometryFile
) const
{
    forAll(polys, i)
    {
        writeEnsDataBinary
        (
            cellFaces[polys[i]].size(),
            ensightGeometryFile
        );
    }
}


void Foam::ensightMesh::writePolysNPointsPerFaceBinary
(
    const labelList& polys,
    const cellList& cellFaces,
    const faceList& faces,
    std::ofstream& ensightGeometryFile
) const
{
    forAll(polys, i)
    {
        const labelList& cf = cellFaces[polys[i]];

        forAll(cf, faceI)
        {
            writeEnsDataBinary
            (
                faces[cf[faceI]].size(),
                ensightGeometryFile
            );
        }
    }
}


void Foam::ensightMesh::writePolysPointsBinary
(
    const labelList& polys,
    const cellList& cellFaces,
    const faceList& faces,
    std::ofstream& ensightGeometryFile
) const
{
    forAll(polys, i)
    {
        const labelList& cf = cellFaces[polys[i]];

        forAll(cf, faceI)
        {
            const face& f = faces[cf[faceI]];

            forAll(f, pointI)
            {
                writeEnsDataBinary(f[pointI] + 1,ensightGeometryFile);
            }
        }
    }
}


void Foam::ensightMesh::writeAllPolysBinary
(
    const labelList& pointToGlobal,
    std::ofstream& ensightGeometryFile
) const
{
    if (meshCellSets_.nPolys)
    {
        const cellList& cellFaces = mesh_.cells();
        // Renumber faces to use global point numbers
        faceList faces(mesh_.faces());
        forAll(faces, i)
        {
            inplaceRenumber(pointToGlobal, faces[i]);
        }

        if (Pstream::master())
        {
            writeEnsDataBinary("nfaced",ensightGeometryFile);
            writeEnsDataBinary(meshCellSets_.nPolys,ensightGeometryFile);
        }

        // Number of faces for each poly cell
        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                toMaster<< meshCellSets_.polys << cellFaces;
            }

            pBufs.finishedSends();

            if (Pstream::master())
            {
                // Master
                writePolysNFacesBinary
                (
                    meshCellSets_.polys,
                    cellFaces,
                    ensightGeometryFile
                );
                // Slaves
                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    UIPstream fromSlave(slave, pBufs);
                    labelList polys(fromSlave);
                    cellList cellFaces(fromSlave);

                    writePolysNFacesBinary
                    (
                        polys,
                        cellFaces,
                        ensightGeometryFile
                    );
                }
            }
        }

        // Number of points for each face of the above list
        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                toMaster<< meshCellSets_.polys << cellFaces << faces;
            }

            pBufs.finishedSends();

            if (Pstream::master())
            {
                // Master
                writePolysNPointsPerFaceBinary
                (
                    meshCellSets_.polys,
                    cellFaces,
                    faces,
                    ensightGeometryFile
                );
                // Slaves
                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    UIPstream fromSlave(slave, pBufs);
                    labelList polys(fromSlave);
                    cellList cellFaces(fromSlave);
                    faceList faces(fromSlave);

                    writePolysNPointsPerFaceBinary
                    (
                        polys,
                        cellFaces,
                        faces,
                        ensightGeometryFile
                    );
                }
            }
        }

        // List of points id for each face of the above list
        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                toMaster<< meshCellSets_.polys << cellFaces << faces;
            }

            pBufs.finishedSends();

            if (Pstream::master())
            {
                // Master
                writePolysPointsBinary
                (
                    meshCellSets_.polys,
                    cellFaces,
                    faces,
                    ensightGeometryFile
                );
                // Slaves
                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    UIPstream fromSlave(slave, pBufs);
                    labelList polys(fromSlave);
                    cellList cellFaces(fromSlave);
                    faceList faces(fromSlave);

                    writePolysPointsBinary
                    (
                        polys,
                        cellFaces,
                        faces,
                        ensightGeometryFile
                    );
                }
            }
        }
    }
}


void Foam::ensightMesh::writeAllPrims
(
    const char* key,
    const label nPrims,
    const cellShapeList& cellShapes,
    OFstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        PstreamBuffers pBufs(Pstream::nonBlocking);

        if (!Pstream::master())
        {
            UOPstream toMaster(Pstream::masterNo(), pBufs);
            toMaster<< cellShapes;
        }

        pBufs.finishedSends();

        if (Pstream::master())
        {
            ensightGeometryFile << key << nl << setw(10) << nPrims << nl;

            writePrims(cellShapes, ensightGeometryFile);

            for (int slave=1; slave<Pstream::nProcs(); slave++)
            {
                UIPstream fromSlave(slave, pBufs);
                cellShapeList cellShapes(fromSlave);

                writePrims(cellShapes, ensightGeometryFile);
            }
        }
    }
}


void Foam::ensightMesh::writeAllPrimsBinary
(
    const char* key,
    const label nPrims,
    const cellShapeList& cellShapes,
    std::ofstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        PstreamBuffers pBufs(Pstream::nonBlocking);

        if (!Pstream::master())
        {
            UOPstream toMaster(Pstream::masterNo(), pBufs);
            toMaster<< cellShapes;
        }

        pBufs.finishedSends();

        if (Pstream::master())
        {
            writeEnsDataBinary(key,ensightGeometryFile);
            writeEnsDataBinary(nPrims,ensightGeometryFile);

            writePrimsBinary(cellShapes, ensightGeometryFile);

            for (int slave=1; slave<Pstream::nProcs(); slave++)
            {
                UIPstream fromSlave(slave, pBufs);
                cellShapeList cellShapes(fromSlave);

                writePrimsBinary(cellShapes, ensightGeometryFile);
            }
        }
    }
}


void Foam::ensightMesh::writeFacePrims
(
    const faceList& patchFaces,
    OFstream& ensightGeometryFile
) const
{
    forAll(patchFaces, i)
    {
        const face& patchFace = patchFaces[i];

        forAll(patchFace, pointI)
        {
            ensightGeometryFile << setw(10) << patchFace[pointI] + 1;
        }
        ensightGeometryFile << nl;
    }
}


void Foam::ensightMesh::writeFacePrimsBinary
(
    const faceList& patchFaces,
    std::ofstream& ensightGeometryFile
) const
{
    forAll(patchFaces, i)
    {
        const face& patchFace = patchFaces[i];

        forAll(patchFace, pointI)
        {
            writeEnsDataBinary(patchFace[pointI] + 1, ensightGeometryFile);
        }
    }
}


void Foam::ensightMesh::writeAllFacePrims
(
    const char* key,
    const labelList& prims,
    const label nPrims,
    const faceList& patchFaces,
    OFstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        PstreamBuffers pBufs(Pstream::nonBlocking);

        if (!Pstream::master())
        {
            UOPstream toMaster(Pstream::masterNo(), pBufs);
            toMaster<< UIndirectList<face>(patchFaces, prims);
        }

        pBufs.finishedSends();

        if (Pstream::master())
        {
            ensightGeometryFile << key << nl << setw(10) << nPrims << nl;

            writeFacePrims
            (
                UIndirectList<face>(patchFaces, prims)(),
                ensightGeometryFile
            );

            for (int slave=1; slave<Pstream::nProcs(); slave++)
            {
                UIPstream fromSlave(slave, pBufs);
                faceList patchFaces(fromSlave);

                writeFacePrims(patchFaces, ensightGeometryFile);
            }
        }
    }
}


void Foam::ensightMesh::writeNSidedNPointsPerFace
(
    const faceList& patchFaces,
    OFstream& ensightGeometryFile
) const
{
    forAll(patchFaces, i)
    {
        ensightGeometryFile << setw(10) << patchFaces[i].size() << nl;
    }
}


void Foam::ensightMesh::writeNSidedPoints
(
    const faceList& patchFaces,
    OFstream& ensightGeometryFile
) const
{
    writeFacePrims(patchFaces, ensightGeometryFile);
}


void Foam::ensightMesh::writeAllNSided
(
    const labelList& prims,
    const label nPrims,
    const faceList& patchFaces,
    OFstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            ensightGeometryFile
                << "nsided" << nl << setw(10) << nPrims << nl;
        }

        // Number of points for each face
        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                toMaster<< UIndirectList<face>(patchFaces, prims);
            }

            pBufs.finishedSends();

            if (Pstream::master())
            {
                writeNSidedNPointsPerFace
                (
                    UIndirectList<face>(patchFaces, prims)(),
                    ensightGeometryFile
                );

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    UIPstream fromSlave(slave, pBufs);
                    faceList patchFaces(fromSlave);

                    writeNSidedNPointsPerFace
                    (
                        patchFaces,
                        ensightGeometryFile
                    );
                }
            }
        }

        // List of points id for each face
        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                toMaster<< UIndirectList<face>(patchFaces, prims);
            }

            pBufs.finishedSends();

            if (Pstream::master())
            {
                writeNSidedPoints
                (
                    UIndirectList<face>(patchFaces, prims)(),
                    ensightGeometryFile
                );

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    UIPstream fromSlave(slave, pBufs);
                    faceList patchFaces(fromSlave);

                    writeNSidedPoints(patchFaces, ensightGeometryFile);
                }
            }
        }
    }
}


void Foam::ensightMesh::writeNSidedPointsBinary
(
    const faceList& patchFaces,
    std::ofstream& ensightGeometryFile
) const
{
    writeFacePrimsBinary(patchFaces, ensightGeometryFile);
}


void Foam::ensightMesh::writeNSidedNPointsPerFaceBinary
(
    const faceList& patchFaces,
    std::ofstream& ensightGeometryFile
) const
{
    forAll(patchFaces, i)
    {
        writeEnsDataBinary(patchFaces[i].size(), ensightGeometryFile);
    }
}


void Foam::ensightMesh::writeAllNSidedBinary
(
    const labelList& prims,
    const label nPrims,
    const faceList& patchFaces,
    std::ofstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            writeEnsDataBinary("nsided",ensightGeometryFile);
            writeEnsDataBinary(nPrims,ensightGeometryFile);
        }

        // Number of points for each face
        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                toMaster<< UIndirectList<face>(patchFaces, prims);
            }

            pBufs.finishedSends();

            if (Pstream::master())
            {
                writeNSidedNPointsPerFaceBinary
                (
                    UIndirectList<face>(patchFaces, prims)(),
                    ensightGeometryFile
                );

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    UIPstream fromSlave(slave, pBufs);
                    faceList patchFaces(fromSlave);

                    writeNSidedNPointsPerFaceBinary
                    (
                        patchFaces,
                        ensightGeometryFile
                    );
                }
            }
        }

        // List of points id for each face
        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                toMaster<< UIndirectList<face>(patchFaces, prims);
            }

            pBufs.finishedSends();

            if (Pstream::master())
            {
                writeNSidedPointsBinary
                (
                    UIndirectList<face>(patchFaces, prims)(),
                    ensightGeometryFile
                );

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    UIPstream fromSlave(slave, pBufs);
                    faceList patchFaces(fromSlave);

                    writeNSidedPointsBinary
                    (
                        patchFaces,
                        ensightGeometryFile
                    );
                }
            }
        }
    }
}


void Foam::ensightMesh::writeAllFacePrimsBinary
(
    const char* key,
    const labelList& prims,
    const label nPrims,
    const faceList& patchFaces,
    std::ofstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        PstreamBuffers pBufs(Pstream::nonBlocking);

        if (!Pstream::master())
        {
            UOPstream toMaster(Pstream::masterNo(), pBufs);
            toMaster<< UIndirectList<face>(patchFaces, prims);
        }

        pBufs.finishedSends();

        if (Pstream::master())
        {
            writeEnsDataBinary(key,ensightGeometryFile);
            writeEnsDataBinary(nPrims,ensightGeometryFile);

            writeFacePrimsBinary
            (
                UIndirectList<face>(patchFaces, prims)(),
                ensightGeometryFile
            );

            for (int slave=1; slave<Pstream::nProcs(); slave++)
            {
                UIPstream fromSlave(slave, pBufs);
                faceList patchFaces(fromSlave);

                writeFacePrimsBinary(patchFaces, ensightGeometryFile);
            }
        }
    }
}


void Foam::ensightMesh::write
(
    const fileName& postProcPath,
    const word& prepend,
    const label timeIndex,
    Ostream& ensightCaseFile
) const
{
    // Find global point numbering
    labelList pointToGlobal;
    labelList uniquePointMap;
    autoPtr<globalIndex> globalPoints = mesh_.globalData().mergePoints
    (
        pointToGlobal,
        uniquePointMap
    );

    const pointField uniquePoints(mesh_.points(), uniquePointMap);

    if (binary_)
    {
        writeBinary
        (
            postProcPath,
            prepend,
            timeIndex,
            ensightCaseFile,
            pointToGlobal,
            uniquePoints,
            globalPoints()
        );
    }
    else
    {
        writeAscii
        (
            postProcPath,
            prepend,
            timeIndex,
            ensightCaseFile,
            pointToGlobal,
            uniquePoints,
            globalPoints()
        );
    }
}


void Foam::ensightMesh::writeAscii
(
    const fileName& postProcPath,
    const word& prepend,
    const label timeIndex,
    Ostream& ensightCaseFile,
    const labelList& pointToGlobal,
    const pointField& uniquePoints,
    const globalIndex& globalPoints
) const
{
    const Time& runTime = mesh_.time();
    //const pointField& points = mesh_.points();
    const cellShapeList& cellShapes = mesh_.cellShapes();



    word timeFile = prepend;

    if (timeIndex == 0)
    {
        timeFile += "000.";
    }
    else if (mesh_.moving())
    {
        timeFile += itoa(timeIndex) + '.';
    }

    // set the filename of the ensight file
    fileName ensightGeometryFileName = timeFile + "mesh";

    OFstream *ensightGeometryFilePtr = NULL;
    if (Pstream::master())
    {
        ensightGeometryFilePtr = new OFstream
        (
            postProcPath/ensightGeometryFileName,
            runTime.writeFormat(),
            runTime.writeVersion(),
            IOstream::UNCOMPRESSED
        );
    }

    OFstream& ensightGeometryFile = *ensightGeometryFilePtr;

    if (Pstream::master())
    {
        // Set Format
        ensightGeometryFile.setf
        (
            ios_base::scientific,
            ios_base::floatfield
        );
        ensightGeometryFile.precision(5);

        ensightGeometryFile
            << "EnSight Geometry File" << nl
            << "written from OpenFOAM-" << Foam::FOAMversion << nl
            << "node id assign" << nl
            << "element id assign" << nl;
    }

    if (patchNames_.empty())
    {
        label nPoints = globalPoints.size();

        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                for (direction d=0; d<vector::nComponents; d++)
                {
                    toMaster<< uniquePoints.component(d);
                }
            }

            pBufs.finishedSends();

            if (Pstream::master())
            {
                ensightGeometryFile
                    << "part" << nl
                    << setw(10) << 1 << nl
                    << "internalMesh" << nl
                    << "coordinates" << nl
                    << setw(10) << nPoints
                    << endl;


                PtrList<UIPstream> fromSlaves(Pstream::nProcs());
                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    fromSlaves.set(slave, new UIPstream(slave, pBufs));
                }

                for (direction d=0; d<vector::nComponents; d++)
                {
                    writePoints(uniquePoints.component(d), ensightGeometryFile);

                    for (int slave=1; slave<Pstream::nProcs(); slave++)
                    {
                        scalarField pointsComponent(fromSlaves[slave]);
                        writePoints(pointsComponent, ensightGeometryFile);
                    }
                }
            }
        }


        writeAllPrims
        (
            "hexa8",
            meshCellSets_.nHexesWedges,
            map         // Rewrite cellShapes to global numbering
            (
                cellShapes,
                meshCellSets_.hexes,
                meshCellSets_.wedges,
                pointToGlobal
            ),
            ensightGeometryFile
        );

        writeAllPrims
        (
            "penta6",
            meshCellSets_.nPrisms,
            map(cellShapes, meshCellSets_.prisms, pointToGlobal),
            ensightGeometryFile
        );

        writeAllPrims
        (
            "pyramid5",
            meshCellSets_.nPyrs,
            map(cellShapes, meshCellSets_.pyrs, pointToGlobal),
            ensightGeometryFile
        );

        writeAllPrims
        (
            "tetra4",
            meshCellSets_.nTets,
            map(cellShapes, meshCellSets_.tets, pointToGlobal),
            ensightGeometryFile
        );

        writeAllPolys
        (
            pointToGlobal,
            ensightGeometryFile
        );
    }


    label ensightPatchI = patchPartOffset_;

    forAll(allPatchNames_, patchi)
    {
        const word& patchName = allPatchNames_[patchi];

        if (patchNames_.empty() || patchNames_.found(patchName))
        {
            const nFacePrimitives& nfp = nPatchPrims_[patchName];

            if (nfp.nTris || nfp.nQuads || nfp.nPolys)
            {
                const polyPatch& p = mesh_.boundaryMesh()[patchi];
                const labelList& tris = boundaryFaceSets_[patchi].tris;
                const labelList& quads = boundaryFaceSets_[patchi].quads;
                const labelList& polys = boundaryFaceSets_[patchi].polys;

                // Renumber the patch points/faces into unique points
                labelList pointToGlobal;
                labelList uniqueMeshPointLabels;
                autoPtr<globalIndex> globalPointsPtr =
                mesh_.globalData().mergePoints
                (
                    p.meshPoints(),
                    p.meshPointMap(),
                    pointToGlobal,
                    uniqueMeshPointLabels
                );

                pointField uniquePoints(mesh_.points(), uniqueMeshPointLabels);
                // Renumber the patch faces
                faceList patchFaces(p.localFaces());
                forAll(patchFaces, i)
                {
                    inplaceRenumber(pointToGlobal, patchFaces[i]);
                }


                {
                    PstreamBuffers pBufs(Pstream::nonBlocking);

                    if (!Pstream::master())
                    {
                        UOPstream toMaster(Pstream::masterNo(), pBufs);
                        for (direction d=0; d<vector::nComponents; d++)
                        {
                            toMaster<< uniquePoints.component(d);
                        }
                    }

                    pBufs.finishedSends();

                    if (Pstream::master())
                    {
                        ensightGeometryFile
                            << "part" << nl
                            << setw(10) << ensightPatchI++ << nl
                            << patchName << nl
                            << "coordinates" << nl
                            << setw(10) << globalPointsPtr().size()
                            << endl;

                        PtrList<UIPstream> fromSlaves(Pstream::nProcs());
                        for (int slave=1; slave<Pstream::nProcs(); slave++)
                        {
                            fromSlaves.set(slave, new UIPstream(slave, pBufs));
                        }

                        for (direction d=0; d<vector::nComponents; d++)
                        {
                            writePoints
                            (
                                uniquePoints.component(d),
                                ensightGeometryFile
                            );

                            for (int slave=1; slave<Pstream::nProcs(); slave++)
                            {
                                scalarField patchPointsComponent
                                (
                                    fromSlaves[slave]
                                );

                                writePoints
                                (
                                    patchPointsComponent,
                                    ensightGeometryFile
                                );
                            }
                        }
                    }
                }

                writeAllFacePrims
                (
                    "tria3",
                    tris,
                    nfp.nTris,
                    patchFaces,
                    ensightGeometryFile
                );

                writeAllFacePrims
                (
                    "quad4",
                    quads,
                    nfp.nQuads,
                    patchFaces,
                    ensightGeometryFile
                );

                writeAllNSided
                (
                    polys,
                    nfp.nPolys,
                    patchFaces,
                    ensightGeometryFile
                );
            }
        }
    }

    if (Pstream::master())
    {
        delete ensightGeometryFilePtr;
    }
}


void Foam::ensightMesh::writeBinary
(
    const fileName& postProcPath,
    const word& prepend,
    const label timeIndex,
    Ostream& ensightCaseFile,
    const labelList& pointToGlobal,
    const pointField& uniquePoints,
    const globalIndex& globalPoints
) const
{
    const cellShapeList& cellShapes = mesh_.cellShapes();

    word timeFile = prepend;

    if (timeIndex == 0)
    {
        timeFile += "000.";
    }
    else if (mesh_.moving())
    {
        timeFile += itoa(timeIndex) + '.';
    }

    // set the filename of the ensight file
    fileName ensightGeometryFileName = timeFile + "mesh";

    std::ofstream *ensightGeometryFilePtr = NULL;

    if (Pstream::master())
    {
        ensightGeometryFilePtr = new std::ofstream
        (
            (postProcPath/ensightGeometryFileName).c_str(),
            ios_base::out | ios_base::binary | ios_base::trunc
        );
        // Check on file opened?
    }

    std::ofstream& ensightGeometryFile = *ensightGeometryFilePtr;

    if (Pstream::master())
    {
        writeEnsDataBinary("C binary", ensightGeometryFile);
        writeEnsDataBinary("EnSight Geometry File", ensightGeometryFile);
        writeEnsDataBinary("written from OpenFOAM", ensightGeometryFile);
        writeEnsDataBinary("node id assign", ensightGeometryFile);
        writeEnsDataBinary("element id assign", ensightGeometryFile);
    }


    if (patchNames_.empty())
    {
        label nPoints = globalPoints.size();

        {
            PstreamBuffers pBufs(Pstream::nonBlocking);

            if (!Pstream::master())
            {
                UOPstream toMaster(Pstream::masterNo(), pBufs);
                for (direction d=0; d<vector::nComponents; d++)
                {
                    toMaster<< uniquePoints.component(d);
                }
            }

            pBufs.finishedSends();

            if (Pstream::master())
            {
                writeEnsDataBinary("part",ensightGeometryFile);
                writeEnsDataBinary(1,ensightGeometryFile);
                writeEnsDataBinary("internalMesh",ensightGeometryFile);
                writeEnsDataBinary("coordinates",ensightGeometryFile);
                writeEnsDataBinary(nPoints,ensightGeometryFile);

                PtrList<UIPstream> fromSlaves(Pstream::nProcs());
                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    fromSlaves.set(slave, new UIPstream(slave, pBufs));
                }

                for (direction d=0; d<vector::nComponents; d++)
                {
                    writeEnsDataBinary
                    (
                        uniquePoints.component(d),
                        ensightGeometryFile
                    );

                    for (int slave=1; slave<Pstream::nProcs(); slave++)
                    {
                        scalarField pointsComponent(fromSlaves[slave]);
                        writeEnsDataBinary
                        (
                            pointsComponent,
                            ensightGeometryFile
                        );
                    }
                }
            }
        }

        writeAllPrimsBinary
        (
            "hexa8",
            meshCellSets_.nHexesWedges,
            map         // Rewrite cellShapes to global numbering
            (
                cellShapes,
                meshCellSets_.hexes,
                meshCellSets_.wedges,
                pointToGlobal
            ),
            ensightGeometryFile
        );

        writeAllPrimsBinary
        (
            "penta6",
            meshCellSets_.nPrisms,
            map(cellShapes, meshCellSets_.prisms, pointToGlobal),
            ensightGeometryFile
        );

        writeAllPrimsBinary
        (
            "pyramid5",
            meshCellSets_.nPyrs,
            map(cellShapes, meshCellSets_.pyrs, pointToGlobal),
            ensightGeometryFile
        );

        writeAllPrimsBinary
        (
            "tetra4",
            meshCellSets_.nTets,
            map(cellShapes, meshCellSets_.tets, pointToGlobal),
            ensightGeometryFile
        );

        writeAllPolysBinary
        (
            pointToGlobal,
            ensightGeometryFile
        );
    }

    label ensightPatchI = patchPartOffset_;
    label iCount = 0;

    forAll(allPatchNames_, patchi)
    {
        iCount ++;
        const word& patchName = allPatchNames_[patchi];

        if (patchNames_.empty() || patchNames_.found(patchName))
        {
            const nFacePrimitives& nfp = nPatchPrims_.find(patchName)();

            if (nfp.nTris || nfp.nQuads || nfp.nPolys)
            {
                const polyPatch& p = mesh_.boundaryMesh()[patchi];
                const labelList& tris = boundaryFaceSets_[patchi].tris;
                const labelList& quads = boundaryFaceSets_[patchi].quads;
                const labelList& polys = boundaryFaceSets_[patchi].polys;

                // Renumber the patch points/faces into unique points
                labelList pointToGlobal;
                labelList uniqueMeshPointLabels;
                autoPtr<globalIndex> globalPointsPtr =
                mesh_.globalData().mergePoints
                (
                    p.meshPoints(),
                    p.meshPointMap(),
                    pointToGlobal,
                    uniqueMeshPointLabels
                );
                pointField uniquePoints(mesh_.points(), uniqueMeshPointLabels);
                // Renumber the patch faces
                faceList patchFaces(p.localFaces());
                forAll(patchFaces, i)
                {
                    inplaceRenumber(pointToGlobal, patchFaces[i]);
                }


                {
                    PstreamBuffers pBufs(Pstream::nonBlocking);

                    if (!Pstream::master())
                    {
                        UOPstream toMaster(Pstream::masterNo(), pBufs);
                        for (direction d=0; d<vector::nComponents; d++)
                        {
                            toMaster<< uniquePoints.component(d);
                        }
                    }

                    pBufs.finishedSends();

                    if (Pstream::master())
                    {
                        writeEnsDataBinary("part",ensightGeometryFile);
                        writeEnsDataBinary(ensightPatchI++,ensightGeometryFile);
                        //writeEnsDataBinary
                        //(patchName.c_str(),ensightGeometryFile);
                        writeEnsDataBinary
                        (   
                            patchName.c_str(),
                            ensightGeometryFile
                        );
                        writeEnsDataBinary("coordinates",ensightGeometryFile);
                        writeEnsDataBinary
                        (
                            globalPointsPtr().size(),
                            ensightGeometryFile
                        );

                        PtrList<UIPstream> fromSlaves(Pstream::nProcs());
                        for (int slave=1; slave<Pstream::nProcs(); slave++)
                        {
                            fromSlaves.set(slave, new UIPstream(slave, pBufs));
                        }

                        for (direction d=0; d<vector::nComponents; d++)
                        {
                            //writePointsBinary
                            writeEnsDataBinary
                            (
                                uniquePoints.component(d),
                                ensightGeometryFile
                            );

                            for (int slave=1; slave<Pstream::nProcs(); slave++)
                            {
                                scalarField patchPointsComponent
                                (
                                    fromSlaves[slave]
                                );

                                //writePointsBinary
                                writeEnsDataBinary
                                (
                                    patchPointsComponent,
                                    ensightGeometryFile
                                );
                            }
                        }
                    }
                }

                writeAllFacePrimsBinary
                (
                    "tria3",
                    tris,
                    nfp.nTris,
                    patchFaces,
                    ensightGeometryFile
                );

                writeAllFacePrimsBinary
                (
                    "quad4",
                    quads,
                    nfp.nQuads,
                    patchFaces,
                    ensightGeometryFile
                );

                writeAllNSidedBinary
                (
                    polys,
                    nfp.nPolys,
                    patchFaces,
                    ensightGeometryFile
                );
            }
        }
    }


    if (Pstream::master())
    {
        delete ensightGeometryFilePtr;
    }
}


// ************************************************************************* //
