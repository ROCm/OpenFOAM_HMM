/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "ensightMesh.H"
#include "fvMesh.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "processorPolyPatch.H"
#include "mapDistribute.H"
#include "stringListOps.H"

#include "ensightFile.H"
#include "ensightGeoFile.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::cellShapeList& Foam::ensightMesh::renumberShapes
(
    cellShapeList& shapes,
    const labelUList& pointToGlobal
)
{
    forAll(shapes, i)
    {
        inplaceRenumber(pointToGlobal, shapes[i]);
    }

    return shapes;
}


Foam::cellShapeList Foam::ensightMesh::map
(
    const cellShapeList& shapes,
    const labelUList& addr,
    const labelUList& pointToGlobal
)
{
    cellShapeList lst(addr.size());

    forAll(addr, i)
    {
        lst[i] = shapes[addr[i]];
        inplaceRenumber(pointToGlobal, lst[i]);
    }

    return lst;
}


void Foam::ensightMesh::writeFaceList
(
    const faceList& faceLst,
    ensightGeoFile& os
)
{
    forAll(faceLst, i)
    {
        const face& f = faceLst[i];

        forAll(f, fp)
        {
            os.write(f[fp] + 1);
        }

        os.newline();
    }
}


void Foam::ensightMesh::writeFaceList
(
    const UIndirectList<face>& faceLst,
    ensightGeoFile& os
)
{
    forAll(faceLst, i)
    {
        const face& f = faceLst[i];

        forAll(f, fp)
        {
            os.write(f[fp] + 1);
        }

        os.newline();
    }
}


void Foam::ensightMesh::writeFaceSizes
(
    const faceList& faceLst,
    ensightGeoFile& os
)
{
    forAll(faceLst, i)
    {
        const face& f = faceLst[i];

        os.write(f.size());
        os.newline();
    }
}


void Foam::ensightMesh::writeFaceSizes
(
    const UIndirectList<face>& faceLst,
    ensightGeoFile& os
)
{
    forAll(faceLst, i)
    {
        const face& f = faceLst[i];

        os.write(f.size());
        os.newline();
    }
}


void Foam::ensightMesh::writeCellShapes
(
    const cellShapeList& shapes,
    ensightGeoFile& os
)
{
    forAll(shapes, i)
    {
        const cellShape& cellPoints = shapes[i];

        // convert global -> local index
        // (note: Ensight indices start with 1)

        // In ASCII, write one cell per line
        forAll(cellPoints, pointI)
        {
            os.write(cellPoints[pointI] + 1);
        }

        os.newline();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightMesh::writePolysNFaces
(
    const labelList& addr,
    const cellList&  cellFaces,
    ensightGeoFile& os
) const
{
    // write the number of faces per element (1/line in ASCII)
    forAll(addr, i)
    {
        const labelUList& cf = cellFaces[addr[i]];

        os.write(cf.size());
        os.newline();
    }
}


void Foam::ensightMesh::writePolysNPointsPerFace
(
    const labelList& addr,
    const cellList& cellFaces,
    const faceList& faces,
    ensightGeoFile& os
) const
{
    // write the number of points per element face (1/line in ASCII)
    forAll(addr, i)
    {
        const labelUList& cf = cellFaces[addr[i]];

        forAll(cf, facei)
        {
            os.write(faces[cf[facei]].size());
            os.newline();
        }
    }
}


void Foam::ensightMesh::writePolysPoints
(
    const labelList& addr,
    const cellList& cellFaces,
    const faceList& faces,
    const labelList& faceOwner,
    ensightGeoFile& os
) const
{
    forAll(addr, i)
    {
        const label cellId = addr[i];
        const labelUList& cf = cellFaces[cellId];

        forAll(cf, facei)
        {
            const label faceId = cf[facei];
            const face& f = faces[faceId];  // face points (in global points)

            if (faceId < faceOwner.size() && faceOwner[faceId] != cellId)
            {
                // internal face, neighbour
                //
                // as per face::reverseFace(), but without copying

                os.write(f[0] + 1);
                for (label ptI = f.size()-1; ptI > 0; --ptI)
                {
                    os.write(f[ptI] + 1);
                }
            }
            else
            {
                forAll(f, ptI)
                {
                    os.write(f[ptI] + 1);
                }
            }

            os.newline();
        }
    }
}


void Foam::ensightMesh::writePolysConnectivity
(
    const labelList& addr,
    const labelList& pointToGlobal,
    ensightGeoFile& os
) const
{
    const cellList&  cellFaces = mesh_.cells();
    const faceList&  meshFaces = mesh_.faces();
    const labelList& faceOwner = mesh_.faceOwner();

    if (Pstream::master())
    {
        // Number of faces for each poly cell

        // Master
        writePolysNFaces(addr, cellFaces, os);

        // Slaves
        for (int slave=1; slave<Pstream::nProcs(); ++slave)
        {
            IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
            labelList addr(fromSlave);
            cellList  cellFaces(fromSlave);

            writePolysNFaces(addr, cellFaces, os);
        }
    }
    else
    {
        OPstream toMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
        toMaster
            << addr
            << cellFaces;
    }

    // Number of points for each face of the above list
    if (Pstream::master())
    {
        // Master
        writePolysNPointsPerFace
        (
            addr,
            cellFaces,
            meshFaces,
            os
        );
        // Slaves
        for (int slave=1; slave<Pstream::nProcs(); ++slave)
        {
            IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
            labelList addr(fromSlave);
            cellList  cellFaces(fromSlave);
            faceList  meshFaces(fromSlave);

            writePolysNPointsPerFace
            (
                addr,
                cellFaces,
                meshFaces,
                os
            );
        }
    }
    else
    {
        OPstream toMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
        toMaster
            << addr
            << cellFaces
            << meshFaces;
    }


    // Renumber faces to use global point numbers
    faceList faces(mesh_.faces());
    forAll(faces, i)
    {
        inplaceRenumber(pointToGlobal, faces[i]);
    }

    // List of points id for each face of the above list
    if (Pstream::master())
    {
        // Master
        writePolysPoints
        (
            addr,
            cellFaces,
            faces,
            faceOwner,
            os
        );
        // Slaves
        for (int slave=1; slave<Pstream::nProcs(); ++slave)
        {
            IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
            labelList addr(fromSlave);
            cellList  cellFaces(fromSlave);
            faceList  faces(fromSlave);
            labelList faceOwner(fromSlave);

            writePolysPoints
            (
                addr,
                cellFaces,
                faces,
                faceOwner,
                os
            );
        }
    }
    else
    {
        OPstream toMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
        toMaster
            << addr
            << cellFaces
            << faces
            << faceOwner;
    }
}


void Foam::ensightMesh::writeCellConnectivity
(
    const ensightCells::elemType elemType,
    const ensightCells& ensCells,
    const labelList& pointToGlobal,
    ensightGeoFile& os
) const
{
    const label nTotal = ensCells.total(elemType);

    if (nTotal)
    {
        const labelUList& addr = ensCells.cellIds(elemType);

        if (Pstream::master())
        {
            os.writeKeyword(ensightCells::key(elemType));
            os.write(nTotal);
            os.newline();
        }

        if (elemType == ensightCells::NFACED)
        {
            writePolysConnectivity
            (
                addr,
                pointToGlobal,
                os
            );
        }
        else
        {
            const cellShapeList shapes = map
            (
                mesh_.cellShapes(),
                addr,
                pointToGlobal
            );


            if (Pstream::master())
            {
                writeCellShapes(shapes, os);

                for (int slave=1; slave<Pstream::nProcs(); ++slave)
                {
                    IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                    cellShapeList received(fromSlave);

                    writeCellShapes(received, os);
                }
            }
            else
            {
                OPstream toMaster
                (
                    Pstream::commsTypes::scheduled,
                    Pstream::masterNo()
                );

                toMaster
                    << shapes;
            }
        }
    }
}


void Foam::ensightMesh::writeCellConnectivity
(
    const ensightCells& ensCells,
    const labelList& pointToGlobal,
    ensightGeoFile& os
) const
{
    for (label typei=0; typei < ensightCells::nTypes; ++typei)
    {
        const ensightCells::elemType what = ensightCells::elemType(typei);

        writeCellConnectivity(what, ensCells, pointToGlobal, os);
    }
}


void Foam::ensightMesh::writeFaceConnectivity
(
    ensightFaces::elemType elemType,
    const label nTotal,
    const faceList& faces,
    ensightGeoFile& os
) const
{
    if (nTotal)
    {
        if (Pstream::master())
        {
            os.writeKeyword(ensightFaces::key(elemType));
            os.write(nTotal);
            os.newline();
        }

        if (elemType == ensightFaces::NSIDED)
        {
            // Number of points per face

            if (Pstream::master())
            {
                writeFaceSizes(faces, os);

                for (int slave=1; slave<Pstream::nProcs(); ++slave)
                {
                    IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                    faceList received(fromSlave);

                    writeFaceSizes(received, os);
                }
            }
            else
            {
                OPstream toMaster
                (
                    Pstream::commsTypes::scheduled,
                    Pstream::masterNo()
                );

                toMaster
                    << faces;
            }
        }


        // List of points id for each face
        if (Pstream::master())
        {
            writeFaceList(faces, os);

            for (int slave=1; slave<Pstream::nProcs(); ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                faceList received(fromSlave);

                writeFaceList(received, os);
            }
        }
        else
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster
                << faces;
        }
    }
}


void Foam::ensightMesh::writeFaceConnectivity
(
    ensightFaces::elemType elemType,
    const label nTotal,
    const faceList& faceLst,
    const labelList& addr,
    ensightGeoFile& os
) const
{
    if (nTotal)
    {
        if (Pstream::master())
        {
            os.writeKeyword(ensightFaces::key(elemType));
            os.write(nTotal);
            os.newline();
        }

        const UIndirectList<face> faces(faceLst, addr);

        if (elemType == ensightFaces::NSIDED)
        {
            // Number of points per face

            if (Pstream::master())
            {
                writeFaceSizes(faces, os);

                for (int slave=1; slave<Pstream::nProcs(); ++slave)
                {
                    IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                    faceList received(fromSlave);

                    writeFaceSizes(received, os);
                }
            }
            else
            {
                OPstream toMaster
                (
                    Pstream::commsTypes::scheduled,
                    Pstream::masterNo()
                );

                toMaster
                    << faces;
            }
        }

        // List of points id per face
        if (Pstream::master())
        {
            writeFaceList(faces, os);

            for (int slave=1; slave<Pstream::nProcs(); ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                faceList received(fromSlave);

                writeFaceList(received, os);
            }
        }
        else
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster
                << faces;
        }
    }
}


void Foam::ensightMesh::writeFaceConnectivity
(
    const ensightFaces& ensFaces,
    const faceList& faceLst,
    ensightGeoFile& os,
    const bool raw
) const
{
    if (raw)
    {
        for (label typei=0; typei < ensightFaces::nTypes; ++typei)
        {
            const ensightFaces::elemType what = ensightFaces::elemType(typei);

            writeFaceConnectivity
            (
                what,
                ensFaces.total(what),
                SubList<face>
                (
                    faceLst,
                    ensFaces.faceIds(what).size(),
                    ensFaces.offset(what)
                ),
                os
            );
        }
    }
    else
    {
        for (label typei=0; typei < ensightFaces::nTypes; ++typei)
        {
            const ensightFaces::elemType what = ensightFaces::elemType(typei);

            writeFaceConnectivity
            (
                what,
                ensFaces.total(what),
                faceLst,
                ensFaces.faceIds(what),
                os
            );
        }
    }
}


void Foam::ensightMesh::writeAllPoints
(
    const label partId,
    const word& ensightPartName,
    const label nPoints,
    const pointField& uniquePoints,
    ensightGeoFile& os
) const
{
    if (Pstream::master())
    {
        os.beginPart(partId, ensightPartName);

        // write points
        os.beginCoordinates(nPoints);

        for (direction cmpt=0; cmpt < point::nComponents; ++cmpt)
        {
            os.writeList(uniquePoints.component(cmpt));

            for (int slave=1; slave<Pstream::nProcs(); ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                scalarField received(fromSlave);
                os.writeList(received);
            }
        }
    }
    else
    {
        for (direction cmpt=0; cmpt < point::nComponents; ++cmpt)
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster
                << uniquePoints.component(cmpt);
        }
    }
}


// ************************************************************************* //
