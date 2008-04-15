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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "ensightMesh.H"
#include "fvMesh.H"
#include "PstreamCombineReduceOps.H"
#include "processorPolyPatch.H"
#include "cellModeller.H"
#include "IOmanip.H"
#include "itoa.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{

class concatPatchNames
{

public:

    void operator()
    (
        HashTable<labelList>& x,
        const HashTable<labelList>& y
    ) const
    {
        for
        (
            HashTable<labelList>::const_iterator iter = y.begin();
            iter != y.end();
            ++iter
        )
        {
            HashTable<labelList>::iterator xiter = x.find(iter.key());

            if (xiter == x.end())
            {
                x.insert(iter.key(), iter());
            }
            else
            {
                labelList& xPatches = xiter();
                const labelList& yPatches = iter();

                label offset = xPatches.size();
                xPatches.setSize(offset + yPatches.size());

                forAll(yPatches, i)
                {
                    xPatches[i + offset] = yPatches[i];
                }
            }
        }
    }
};

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from fvMesh
Foam::ensightMesh::ensightMesh(const fvMesh& fMesh, const argList& args)
:
    mesh(fMesh),
    meshCellSets(mesh.nCells()),
    boundaryFaceSets(mesh.boundary().size())
{
    forAll (mesh.boundaryMesh(), patchi)
    {
        if (typeid(mesh.boundaryMesh()[patchi]) != typeid(processorPolyPatch))
        {
            if (!allPatchNames.found(mesh.boundaryMesh()[patchi].name()))
            {
                allPatchNames.insert
                (
                    mesh.boundaryMesh()[patchi].name(),
                    labelList(1, Pstream::myProcNo())
                );

                patchIndices.insert
                (
                    mesh.boundaryMesh()[patchi].name(),
                    patchi
                );
            }
        }
    }

    combineReduce(allPatchNames, concatPatchNames());

    if (args.options().found("patches"))
    {
        wordList patchNameList(IStringStream(args.options()["patches"])());

        if (!patchNameList.size())
        {
            patchNameList = allPatchNames.toc();
        }

        forAll (patchNameList, i)
        {
            patchNames.insert(patchNameList[i]);
        }
    }


    const cellShapeList& cellShapes = mesh.cellShapes();

    const cellModel& tet = *(cellModeller::lookup("tet"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& wedge = *(cellModeller::lookup("wedge"));
    const cellModel& hex = *(cellModeller::lookup("hex"));

    labelList& tets = meshCellSets.tets;
    labelList& pyrs = meshCellSets.pyrs;
    labelList& prisms = meshCellSets.prisms;
    labelList& wedges = meshCellSets.wedges;
    labelList& hexes = meshCellSets.hexes;
    labelList& polys = meshCellSets.polys;

    // Count the shapes
    label nTets = 0;
    label nPyrs = 0;
    label nPrisms = 0;
    label nWedges = 0;
    label nHexes = 0;
    label nPolys = 0;

    if (!patchNames.size())
    {
        forAll(cellShapes, celli)
        {
            const cellShape& cellShape = cellShapes[celli];
            const cellModel& cellModel = cellShape.model();

            if (cellModel == tet)
            {
                tets[nTets++] = celli;
            }
            else if (cellModel == pyr)
            {
                pyrs[nPyrs++] = celli;
            }
            else if (cellModel == prism)
            {
                prisms[nPrisms++] = celli;
            }
            else if (cellModel == wedge)
            {
                wedges[nWedges++] = celli;
            }
            else if (cellModel == hex)
            {
                hexes[nHexes++] = celli;
            }
            else
            {
                polys[nPolys++] = celli;
            }
        }

        tets.setSize(nTets);
        pyrs.setSize(nPyrs);
        prisms.setSize(nPrisms);
        wedges.setSize(nWedges);
        hexes.setSize(nHexes);
        polys.setSize(nPolys);

        meshCellSets.nTets = nTets;
        reduce(meshCellSets.nTets, sumOp<label>());

        meshCellSets.nPyrs = nPyrs;
        reduce(meshCellSets.nPyrs, sumOp<label>());

        meshCellSets.nPrisms = nPrisms;
        reduce(meshCellSets.nPrisms, sumOp<label>());

        meshCellSets.nHexesWedges = nHexes + nWedges;
        reduce(meshCellSets.nHexesWedges, sumOp<label>());

        meshCellSets.nPolys = nPolys;
        reduce(meshCellSets.nPolys, sumOp<label>());
    }


    forAll (mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].size())
        {
            const polyPatch& p = mesh.boundaryMesh()[patchi];

            labelList& tris = boundaryFaceSets[patchi].tris;
            labelList& quads = boundaryFaceSets[patchi].quads;
            labelList& polys = boundaryFaceSets[patchi].polys;

            tris.setSize(p.size());
            quads.setSize(p.size());
            polys.setSize(p.size());

            label nTris = 0;
            label nQuads = 0;
            label nPolys = 0;

            forAll(p, facei)
            {
                const face& f = p[facei];

                if (f.size() == 3)
                {
                    tris[nTris++] = facei;
                }
                else if (f.size() == 4)
                {
                    quads[nQuads++] = facei;
                }
                else
                {
                    polys[nPolys++] = facei;
                }
            }

            tris.setSize(nTris);
            quads.setSize(nQuads);
            polys.setSize(nPolys);
        }
    }

    for
    (
        HashTable<labelList>::const_iterator iter = allPatchNames.begin();
        iter != allPatchNames.end();
        ++iter
    )
    {
        const word& patchName = iter.key();
        nFacePrims nfp;

        if (!patchNames.size() || patchNames.found(patchName))
        {
            if (patchIndices.found(patchName))
            {
                label patchi = patchIndices.find(patchName)();

                nfp.nPoints = mesh.boundaryMesh()[patchi].localPoints().size();
                nfp.nTris = boundaryFaceSets[patchi].tris.size();
                nfp.nQuads = boundaryFaceSets[patchi].quads.size();
                nfp.nPolys = boundaryFaceSets[patchi].polys.size();
            }
        }

        reduce(nfp.nPoints, sumOp<label>());
        reduce(nfp.nTris, sumOp<label>());
        reduce(nfp.nQuads, sumOp<label>());
        reduce(nfp.nPolys, sumOp<label>());

        nPatchPrims.insert(patchName, nfp);
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
    forAll(pointsComponent, pointi)
    {
        ensightGeometryFile << setw(12) << float(pointsComponent[pointi]) << nl;
    }
}


Foam::cellShapeList Foam::ensightMesh::map
(
    const cellShapeList& cellShapes,
    const labelList& prims
) const
{
    cellShapeList mcsl(prims.size());

    forAll(prims, i)
    {
        mcsl[i] = cellShapes[prims[i]];
    }

    return mcsl;
}


Foam::cellShapeList Foam::ensightMesh::map
(
    const cellShapeList& cellShapes,
    const labelList& hexes,
    const labelList& wedges
) const
{
    cellShapeList mcsl(hexes.size() + wedges.size());

    forAll(hexes, i)
    {
        mcsl[i] = cellShapes[hexes[i]];
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
    }

    return mcsl;
}


void Foam::ensightMesh::writePrims
(
    const cellShapeList& cellShapes,
    const label pointOffset,
    OFstream& ensightGeometryFile
) const
{
    label po = pointOffset + 1;

    forAll(cellShapes, i)
    {
        const cellShape& cellPoints = cellShapes[i];

        forAll(cellPoints, pointi)
        {
            ensightGeometryFile<< setw(10) << cellPoints[pointi] + po;
        }
        ensightGeometryFile << nl;
    }
}


void Foam::ensightMesh::writePolys
(
    const labelList& polys,
    const cellList& cellFaces,
    const faceList& faces,
    const label pointOffset,
    OFstream& ensightGeometryFile
) const
{
    if (polys.size())
    {
        ensightGeometryFile
            << "nfaced" << nl << setw(10) << polys.size() << nl;

        label po = pointOffset + 1;

        forAll(polys, i)
        {
            ensightGeometryFile
                << setw(10) << cellFaces[polys[i]].size() << nl;
        }

        forAll(polys, i)
        {
            const labelList& cf = cellFaces[polys[i]];

            forAll(cf, facei)
            {
                ensightGeometryFile
                    << setw(10) << faces[cf[facei]].size() << nl;
            }
        }

        forAll(polys, i)
        {
            const labelList& cf = cellFaces[polys[i]];

            forAll(cf, facei)
            {
                const face& f = faces[cf[facei]];

                forAll(f, pointi)
                {
                    ensightGeometryFile << setw(10) << f[pointi] + po;
                }
                ensightGeometryFile << nl;
            }
        }
    }
}


void Foam::ensightMesh::writeAllPrims
(
    const char* key,
    const label nPrims,
    const cellShapeList& cellShapes,
    const labelList& pointOffsets,
    OFstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            ensightGeometryFile << key << nl << setw(10) << nPrims << nl;

            writePrims(cellShapes, 0, ensightGeometryFile);

            for (int slave=1; slave<Pstream::nProcs(); slave++)
            {
                IPstream fromSlave(Pstream::scheduled, slave);
                cellShapeList cellShapes(fromSlave);

                writePrims
                (
                    cellShapes,
                    pointOffsets[slave-1],
                    ensightGeometryFile
                );
            }
        }
        else
        {
            OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
            toMaster<< cellShapes;
        }
    }
}


void Foam::ensightMesh::writeFacePrims
(
    const char* key,
    const faceList& patchFaces,
    const label pointOffset,
    OFstream& ensightGeometryFile
) const
{
    if (patchFaces.size())
    {
        if (word(key) == "nsided")
        {
            ensightGeometryFile
                << key << nl << setw(10) << patchFaces.size() << nl;

            forAll(patchFaces, i)
            {
                ensightGeometryFile
                    << setw(10) << patchFaces[i].size() << nl;
            }
        }

        label po = pointOffset + 1;

        forAll(patchFaces, i)
        {
            const face& patchFace = patchFaces[i];

            forAll(patchFace, pointi)
            {
                ensightGeometryFile << setw(10) << patchFace[pointi] + po;
            }
            ensightGeometryFile << nl;
        }
    }
}


Foam::faceList Foam::ensightMesh::map
(
    const faceList& patchFaces,
    const labelList& prims
) const
{
    faceList ppf(prims.size());

    forAll (prims, i)
    {
        ppf[i] = patchFaces[prims[i]];
    }

    return ppf;
}


void Foam::ensightMesh::writeAllFacePrims
(
    const char* key,
    const labelList& prims,
    const label nPrims,
    const faceList& patchFaces,
    const labelList& pointOffsets,
    const labelList& patchProcessors,
    OFstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            if (word(key) != "nsided")
            {
                ensightGeometryFile << key << nl << setw(10) << nPrims << nl;
            }

            if (&prims != NULL)
            {
                writeFacePrims
                (
                    key,
                    map(patchFaces, prims),
                    0,
                    ensightGeometryFile
                );
            }

            forAll (patchProcessors, i)
            {
                if (patchProcessors[i] != 0)
                {
                    label slave = patchProcessors[i];
                    IPstream fromSlave(Pstream::scheduled, slave);
                    faceList patchFaces(fromSlave);

                    writeFacePrims
                    (
                        key,
                        patchFaces,
                        pointOffsets[i],
                        ensightGeometryFile
                    );
                }
            }
        }
        else if (&prims != NULL)
        {
            OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
            toMaster<< map(patchFaces, prims);
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
    const Time& runTime = mesh.time();
    const pointField& points = mesh.points();
    const cellList& cellFaces = mesh.cells();
    const faceList& faces = mesh.faces();
    const cellShapeList& cellShapes = mesh.cellShapes();

    word timeFile = prepend;

    if (timeIndex == 0)
    {
        timeFile += "000.";
    }
    else if (mesh.moving())
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
            runTime.writeCompression()
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
            << "OpenFOAM Geometry File " << nl
            << "EnSight 8.2.6" << nl
            << "node id assign" << nl
            << "element id assign" << nl;
    }

    labelList pointOffsets(Pstream::nProcs(), 0);

    if (!patchNames.size())
    {
        label nPoints = points.size();
        Pstream::gather(nPoints, sumOp<label>());

        if (Pstream::master())
        {
            ensightGeometryFile
                << "part" << nl
                << setw(10) << 1 << nl
                << "FOAM cells" << nl
                << "coordinates" << nl
                << setw(10) << nPoints
                << endl;

            for (direction d=0; d<vector::nComponents; d++)
            {
                writePoints(points.component(d), ensightGeometryFile);
                pointOffsets[0] = points.size();

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    scalarField pointsComponent(fromSlave);
                    writePoints(pointsComponent, ensightGeometryFile);
                    pointOffsets[slave] =
                        pointOffsets[slave-1]
                      + pointsComponent.size();
                }
            }
        }
        else
        {
            for (direction d=0; d<vector::nComponents; d++)
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< points.component(d);
            }
        }

        writeAllPrims
        (
            "hexa8",
            meshCellSets.nHexesWedges,
            map(cellShapes, meshCellSets.hexes, meshCellSets.wedges),
            pointOffsets,
            ensightGeometryFile
        );

        writeAllPrims
        (
            "penta6",
            meshCellSets.nPrisms,
            map(cellShapes, meshCellSets.prisms),
            pointOffsets,
            ensightGeometryFile
        );

        writeAllPrims
        (
            "pyramid5",
            meshCellSets.nPyrs,
            map(cellShapes, meshCellSets.pyrs),
            pointOffsets,
            ensightGeometryFile
        );

        writeAllPrims
        (
            "tetra4",
            meshCellSets.nTets,
            map(cellShapes, meshCellSets.tets),
            pointOffsets,
            ensightGeometryFile
        );


        if (meshCellSets.nPolys)
        {
            if (Pstream::master())
            {
                /*
                ensightGeometryFile
                    << "nfaced" << nl
                    << setw(10) << meshCellSets.nPolys << nl;
                */
                writePolys
                (
                    meshCellSets.polys,
                    cellFaces,
                    faces,
                    0,
                    ensightGeometryFile
                );

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    labelList polys(fromSlave);
                    cellList cellFaces(fromSlave);
                    faceList faces(fromSlave);

                    writePolys
                    (
                        polys,
                        cellFaces,
                        faces,
                        pointOffsets[slave-1],
                        ensightGeometryFile
                    );
                }
            }
            else
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< meshCellSets.polys << cellFaces << faces;
            }
        }
    }


    label ensightPatchi = 2;

    for
    (
        HashTable<labelList>::const_iterator iter = allPatchNames.begin();
        iter != allPatchNames.end();
        ++iter
    )
    {
        const labelList& patchProcessors = iter();

        if (!patchNames.size() || patchNames.found(iter.key()))
        {
            const word& patchName = iter.key();
            const nFacePrims& nfp = nPatchPrims.find(patchName)();

            const labelList *trisPtr = NULL;
            const labelList *quadsPtr = NULL;
            const labelList *polysPtr = NULL;

            const pointField *patchPointsPtr = NULL;
            const faceList *patchFacesPtr = NULL;

            if (patchIndices.found(iter.key()))
            {
                label patchi = patchIndices.find(iter.key())();
                const polyPatch& p = mesh.boundaryMesh()[patchi];

                trisPtr = &boundaryFaceSets[patchi].tris;
                quadsPtr = &boundaryFaceSets[patchi].quads;
                polysPtr = &boundaryFaceSets[patchi].polys;

                patchPointsPtr = &(p.localPoints());
                patchFacesPtr = &(p.localFaces());
            }

            const labelList& tris = *trisPtr;
            const labelList& quads = *quadsPtr;
            const labelList& polys = *polysPtr;
            const pointField& patchPoints = *patchPointsPtr;
            const faceList& patchFaces = *patchFacesPtr;

            if (nfp.nTris || nfp.nQuads || nfp.nPolys)
            {
                labelList patchPointOffsets(Pstream::nProcs(), 0);

                if (Pstream::master())
                {
                    ensightGeometryFile
                        << "part" << nl
                        << setw(10) << ensightPatchi++ << nl
                        << patchName << nl
                        << "coordinates" << nl
                        << setw(10) << nfp.nPoints
                        << endl;

                    for (direction d=0; d<vector::nComponents; d++)
                    {
                        if (patchPointsPtr)
                        {
                            writePoints
                            (
                                patchPoints.component(d),
                                ensightGeometryFile
                            );
                        }

                        patchPointOffsets = 0;

                        forAll (patchProcessors, i)
                        {
                            if (patchProcessors[i] != 0)
                            {
                                label slave = patchProcessors[i];
                                IPstream fromSlave(Pstream::scheduled, slave);
                                scalarField patchPointsComponent(fromSlave);

                                writePoints
                                (
                                    patchPointsComponent,
                                    ensightGeometryFile
                                );

                                if (i < Pstream::nProcs()-1)
                                {
                                    patchPointOffsets[i+1] =
                                        patchPointOffsets[i]
                                      + patchPointsComponent.size();
                                }
                            }
                            else
                            {
                                if (i < Pstream::nProcs()-1)
                                {
                                    patchPointOffsets[i+1] =
                                        patchPointOffsets[i]
                                      + patchPoints.size();
                                }
                            }
                        }
                    }
                }
                else if (patchPointsPtr)
                {
                    for (direction d=0; d<vector::nComponents; d++)
                    {
                        OPstream toMaster
                        (
                            Pstream::scheduled,
                            Pstream::masterNo()
                        );
                        toMaster<< patchPoints.component(d);
                    }
                }

                writeAllFacePrims
                (
                    "tria3",
                    tris,
                    nfp.nTris,
                    patchFaces,
                    patchPointOffsets,
                    patchProcessors,
                    ensightGeometryFile
                );

                writeAllFacePrims
                (
                    "quad4",
                    quads,
                    nfp.nQuads,
                    patchFaces,
                    patchPointOffsets,
                    patchProcessors,
                    ensightGeometryFile
                );

                writeAllFacePrims
                (
                    "nsided",
                    polys,
                    nfp.nPolys,
                    patchFaces,
                    patchPointOffsets,
                    patchProcessors,
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
