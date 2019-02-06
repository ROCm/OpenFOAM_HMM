/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    for (cellShape& shape : shapes)
    {
        inplaceRenumber(pointToGlobal, shape);
    }

    return shapes;
}


Foam::cellShapeList Foam::ensightMesh::renumberShapes
(
    const cellShapeList& shapes,
    const labelUList& addr,
    const labelUList& pointToGlobal
)
{
    cellShapeList list(shapes, addr);

    renumberShapes(list, pointToGlobal);

    return list;
}


void Foam::ensightMesh::writeFaceList
(
    const faceList& faceLst,
    ensightGeoFile& os
)
{
    for (const face& f : faceLst)
    {
        for (const label labi : f)
        {
            os.write(labi + 1);
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
    for (const face& f : faceLst)
    {
        for (const label labi : f)
        {
            os.write(labi + 1);
        }

        os.newline();
    }
}


Foam::labelList Foam::ensightMesh::getFaceSizes
(
    const faceList& faceLst
)
{
    labelList list(faceLst.size());

    auto outIter = list.begin();

    for (const face& f : faceLst)
    {
        *outIter = f.size();
        ++outIter;
    }

    return list;
}


Foam::labelList Foam::ensightMesh::getFaceSizes
(
    const UIndirectList<face>& faceLst
)
{
    labelList list(faceLst.size());

    auto outIter = list.begin();

    for (const face& f : faceLst)
    {
        *outIter = f.size();
        ++outIter;
    }

    return list;
}


void Foam::ensightMesh::writeFaceSizes
(
    const faceList& faceLst,
    ensightGeoFile& os
)
{
    for (const face& f : faceLst)
    {
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
    for (const face& f : faceLst)
    {
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
    for (const cellShape& cellPoints : shapes)
    {
        // Convert global -> local index
        // (note: Ensight indices start with 1)

        // In ASCII, write one cell per line
        for (const label pointi : cellPoints)
        {
            os.write(pointi + 1);
        }

        os.newline();
    }
}


Foam::labelList Foam::ensightMesh::getPolysNFaces
(
    const labelUList& addr,
    const cellList& cellFaces
)
{
    labelList list(addr.size());

    auto outIter = list.begin();

    // The number of faces per element
    for (const label cellId : addr)
    {
        const labelUList& cf = cellFaces[cellId];

        *outIter = cf.size();
        ++outIter;
    }

    return list;
}


void Foam::ensightMesh::writePolysNFaces
(
    const labelUList& addr,
    const cellList& cellFaces,
    ensightGeoFile& os
)
{
    // Write the number of faces per element (1/line in ASCII)
    for (const label cellId : addr)
    {
        const labelUList& cf = cellFaces[cellId];

        os.write(cf.size());
        os.newline();
    }
}


Foam::labelList Foam::ensightMesh::getPolysNPointsPerFace
(
    const labelUList& addr,
    const cellList& cellFaces,
    const faceList& faces
)
{
    // Count the number of faces per element

    label nTotFaces = 0;
    for (const label cellId : addr)
    {
        const labelUList& cf = cellFaces[cellId];

        nTotFaces += cf.size();
    }

    labelList list(nTotFaces);

    auto outIter = list.begin();

    // The number of points per element face
    for (const label cellId : addr)
    {
        const labelUList& cf = cellFaces[cellId];

        for (const label facei : cf)
        {
            *outIter = faces[facei].size();
            ++outIter;
        }
    }

    return list;
}


void Foam::ensightMesh::writePolysNPointsPerFace
(
    const labelUList& addr,
    const cellList& cellFaces,
    const faceList& faces,
    ensightGeoFile& os
)
{
    // Write the number of points per element face (1/line in ASCII)
    for (const label cellId : addr)
    {
        const labelUList& cf = cellFaces[cellId];

        for (const label facei : cf)
        {
            os.write(faces[facei].size());
            os.newline();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightMesh::writePolysPoints
(
    const labelUList& addr,
    const cellList& cellFaces,
    const faceList& faces,
    const labelList& faceOwner,
    ensightGeoFile& os
)
{
    for (const label cellId : addr)
    {
        const labelUList& cf = cellFaces[cellId];

        for (const label faceId : cf)
        {
            const face& f = faces[faceId];  // face points (in global points)

            if (faceId < faceOwner.size() && faceOwner[faceId] != cellId)
            {
                // internal face, neighbour
                //
                // as per face::reverseFace(), but without copying

                os.write(f[0] + 1);
                for (label pti = f.size()-1; pti > 0; --pti)
                {
                    os.write(f[pti] + 1);
                }
            }
            else
            {
                for (const label labi : f)
                {
                    os.write(labi + 1);
                }
            }

            os.newline();
        }
    }
}


void Foam::ensightMesh::writePolysConnectivity
(
    const labelUList& addr,
    const labelList& pointToGlobal,
    ensightGeoFile& os
) const
{
    const cellList&  cellFaces = mesh_.cells();
    const faceList&  meshFaces = mesh_.faces();
    const labelList& faceOwner = mesh_.faceOwner();

    // Number of faces for each poly cell
    if (Pstream::master())
    {
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
    for (face& f : faces)
    {
        inplaceRenumber(pointToGlobal, f);
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
            const cellShapeList shapes
            (
                renumberShapes
                (
                    mesh_.cellShapes(),
                    addr,
                    pointToGlobal
                )
            );


            if (Pstream::master())
            {
                writeCellShapes(shapes, os);

                for (int slave=1; slave<Pstream::nProcs(); ++slave)
                {
                    IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                    cellShapeList recv(fromSlave);

                    writeCellShapes(recv, os);
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
        const ensightCells::elemType what =
            ensightCells::elemType(typei);

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
                    faceList recv(fromSlave);

                    writeFaceSizes(recv, os);
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
                faceList recv(fromSlave);

                writeFaceList(recv, os);
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
    const labelUList& addr,
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
                    faceList recv(fromSlave);

                    writeFaceSizes(recv, os);
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
                faceList recv(fromSlave);

                writeFaceList(recv, os);
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
    for (label typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const ensightFaces::elemType what =
            ensightFaces::elemType(typei);

        if (raw)
        {
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
        else
        {
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

        // Write points
        os.beginCoordinates(nPoints);

        for (direction cmpt=0; cmpt < point::nComponents; ++cmpt)
        {
            os.writeList(uniquePoints.component(cmpt));

            for (int slave=1; slave<Pstream::nProcs(); ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                scalarField recv(fromSlave);
                os.writeList(recv);
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
