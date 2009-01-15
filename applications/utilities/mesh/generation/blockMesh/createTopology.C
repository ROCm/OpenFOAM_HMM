/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "blockMesh.H"
#include "Time.H"
#include "preservePatchTypes.H"
#include "emptyPolyPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::blockMesh::blockLabelsOK
(
    const label blockLabel,
    const pointField& points,
    const cellShape& blockShape
)
{
    bool ok = true;

    for (int i=0; i<blockShape.size(); i++)
    {
        if (blockShape[i] < 0)
        {
            ok = false;

            WarningIn
            (
                "bool Foam::blockMesh::blockLabelsOK"
                "(const label blockLabel, const pointField& points, "
                "const cellShape& blockShape)"
            )   << "block " << blockLabel
                << " point label " << blockShape[i]
                << " less than zero" << endl;
        }
        else if (blockShape[i] >= points.size())
        {
            ok = false;

            WarningIn
            (
                "bool Foam::blockMesh::blockLabelsOK"
                "(const label blockLabel, const pointField& points, "
                "const cellShape& blockShape)"
            )   << "block " << blockLabel
                << " point label " << blockShape[i]
                << " larger than " << points.size() - 1
                << " the largest defined point label" << endl;
        }
    }

    return ok;
}


bool Foam::blockMesh::patchLabelsOK
(
    const label patchLabel,
    const pointField& points,
    const faceList& patchFaces
)
{
    bool ok = true;

    for (label facei=0; facei<patchFaces.size(); facei++)
    {
        const labelList& labels = patchFaces[facei];

        for (int i=0; i<labels.size(); i++)
        {
            if (labels[i] < 0)
            {
                ok = false;

                WarningIn
                (
                    "bool Foam::blockMesh::patchLabelsOK"
                    "(const label patchLabel, const pointField& points, "
                    "const faceList& patchFaces)"
                )   << "patch " << patchLabel
                    << " face " << facei
                    << " point label " << labels[i]
                    << " less than zero" << endl;
            }
            else if (labels[i] >= points.size())
            {
                ok = false;

                WarningIn
                (
                    "bool Foam::blockMesh::patchLabelsOK"
                    "(const label patchLabel, const pointField& points, "
                    "const faceList& patchFaces)"
                )   << "patch " << patchLabel
                    << " face " << facei
                    << " point label " << labels[i]
                    << " larger than " << points.size() - 1
                    << " the largest defined point label" << endl;
            }
        }
    }

    return ok;
}


bool Foam::blockMesh::readPatches
(
    const dictionary& meshDescription,
    const pointField& tmpBlockPoints,
    faceListList& tmpBlocksPatches,
    wordList& patchNames,
    wordList& patchTypes,
    wordList& nbrPatchNames
)
{
    bool topologyOK = true;

    ITstream& patchStream(meshDescription.lookup("patches"));

    // read number of patches in mesh
    label nPatches = 0;

    token firstToken(patchStream);

    if (firstToken.isLabel())
    {
        nPatches = firstToken.labelToken();

        tmpBlocksPatches.setSize(nPatches);
        patchNames.setSize(nPatches);
        patchTypes.setSize(nPatches);
        nbrPatchNames.setSize(nPatches);
    }
    else
    {
        patchStream.putBack(firstToken);
    }

    // Read beginning of blocks
    patchStream.readBegin("patches");

    nPatches = 0;

    token lastToken(patchStream);
    while
    (
        !(
            lastToken.isPunctuation()
            && lastToken.pToken() == token::END_LIST
        )
    )
    {
        if (tmpBlocksPatches.size() <= nPatches)
        {
            tmpBlocksPatches.setSize(nPatches + 1);
            patchNames.setSize(nPatches + 1);
            patchTypes.setSize(nPatches + 1);
            nbrPatchNames.setSize(nPatches + 1);
        }

        patchStream.putBack(lastToken);

        patchStream
            >> patchTypes[nPatches]
            >> patchNames[nPatches];

        // Read optional neighbour patch name
        if (patchTypes[nPatches] == cyclicPolyPatch::typeName)
        {
            patchStream >> lastToken;
            if (lastToken.isWord())
            {
                nbrPatchNames[nPatches] = lastToken.wordToken();
            }
            else
            {
                patchStream.putBack(lastToken);
            }
        }

        // Read patch faces
        patchStream >> tmpBlocksPatches[nPatches];


        // Catch multiple patches asap.
        for (label i = 0; i < nPatches; i++)
        {
            if (patchNames[nPatches] == patchNames[i])
            {
                FatalErrorIn
                (
                    "blockMesh::createTopology(IOdictionary&)"
                )   << "Duplicate patch " << patchNames[nPatches]
                    << " at line " << patchStream.lineNumber()
                    << ". Exiting !" << nl
                    << exit(FatalError);
            }
        }

        topologyOK = topologyOK && patchLabelsOK
        (
            nPatches,
            tmpBlockPoints,
            tmpBlocksPatches[nPatches]
        );

        nPatches++;


        // Split old style cyclics

        if (patchTypes[nPatches-1] == cyclicPolyPatch::typeName)
        {
            if (nbrPatchNames[nPatches] == word::null)
            {
                word halfA = patchNames[nPatches-1] + "_half0";
                word halfB = patchNames[nPatches-1] + "_half1";

                WarningIn("blockMesh::createTopology(IOdictionary&)")
                    << "Old-style cyclic definition."
                    << " Splitting patch "
                    << patchNames[nPatches-1] << " into two halves "
                    << halfA << " and " << halfB << endl
                    << "    Alternatively use new syntax "
                    << " cyclic <name> <neighbourname> <faces>" << endl;

                // Add extra patch
                if (tmpBlocksPatches.size() <= nPatches)
                {
                    tmpBlocksPatches.setSize(nPatches + 1);
                    patchNames.setSize(nPatches + 1);
                    patchTypes.setSize(nPatches + 1);
                    nbrPatchNames.setSize(nPatches + 1);
                }

                patchNames[nPatches-1] = halfA;
                nbrPatchNames[nPatches-1] = halfB;
                patchTypes[nPatches] = patchTypes[nPatches-1];
                patchNames[nPatches] = halfB;
                nbrPatchNames[nPatches] = halfA;

                // Split faces
                if ((tmpBlocksPatches[nPatches-1].size() % 2) != 0)
                {
                    FatalErrorIn
                    (
                        "blockMesh::createTopology(IOdictionary&)"
                    )   << "Size of cyclic faces is not a multiple of 2 :"
                        << tmpBlocksPatches[nPatches-1]
                        << exit(FatalError);
                }
                label sz = tmpBlocksPatches[nPatches-1].size()/2;
                faceList unsplitFaces(tmpBlocksPatches[nPatches-1], true);
                tmpBlocksPatches[nPatches-1] = faceList
                (
                    SubList<face>(unsplitFaces, sz)
                );
                tmpBlocksPatches[nPatches] = faceList
                (
                    SubList<face>(unsplitFaces, sz, sz)
                );

                nPatches++;
            }
        }

        patchStream >> lastToken;
    }
    patchStream.putBack(lastToken);

    // Read end of blocks
    patchStream.readEnd("patches");

    return topologyOK;
}


bool Foam::blockMesh::readBoundary
(
    const dictionary& meshDescription,
    const pointField& tmpBlockPoints,
    faceListList& tmpBlocksPatches,
    PtrList<dictionary>& patchDicts
)
{
    bool topologyOK = true;

    // Read like boundary file
    const PtrList<entry> patchesInfo
    (
        meshDescription.lookup("boundary")
    );

    tmpBlocksPatches.setSize(patchesInfo.size());
    patchDicts.setSize(patchesInfo.size());

    forAll(tmpBlocksPatches, patchI)
    {
        const entry& patchInfo = patchesInfo[patchI];

        // Construct dictionary and add name
        patchDicts.set(patchI, new dictionary(patchInfo.dict()));
        patchDicts[patchI].set("name", patchInfo.keyword());
        // Read block faces
        patchDicts[patchI].lookup("faces") >> tmpBlocksPatches[patchI];

        topologyOK = topologyOK && patchLabelsOK
        (
            patchI,
            tmpBlockPoints,
            tmpBlocksPatches[patchI]
        );
    }

    return topologyOK;
}


void Foam::blockMesh::createCellShapes
(
    const pointField& tmpBlockPoints,
    PtrList<cellShape>& tmpBlockCells
)
{
    const blockMesh& blocks = *this;

    tmpBlockCells.setSize(blocks.size());
    forAll(blocks, blockLabel)
    {
        tmpBlockCells.set
        (
            blockLabel,
            new cellShape(blocks[blockLabel].blockDef().blockShape())
        );

        if (tmpBlockCells[blockLabel].mag(tmpBlockPoints) < 0.0)
        {
            WarningIn
            (
                "blockMesh::createTopology(IOdictionary&)"
            )   << "negative volume block : " << blockLabel
                << ", probably defined inside-out" << endl;
        }
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::polyMesh* Foam::blockMesh::createTopology(IOdictionary& meshDescription)
{
    bool topologyOK = true;

    blockMesh& blocks = *this;

    word defaultPatchName = "defaultFaces";
    word defaultPatchType = emptyPolyPatch::typeName;

    // get names and types for the unassigned patch faces
    if (meshDescription.found("defaultPatch"))
    {
        const dictionary& defaultPatch =
            meshDescription.subDict("defaultPatch");

        // this is a bit heavy handed (and ugly), but there is currently
        // no easy way to rename polyMesh patches subsequently
        if (defaultPatch.found("name"))
        {
            defaultPatch.lookup("name") >> defaultPatchName;
        }

        if (defaultPatch.found("type"))
        {
            defaultPatch.lookup("type") >> defaultPatchType;
        }
    }

    Info<< nl << "Creating blockCorners" << endl;

    // create blockCorners
    pointField tmpBlockPoints(meshDescription.lookup("vertices"));

    if (meshDescription.found("edges"))
    {
        // read number of non-linear edges in mesh
        Info<< nl << "Creating curved edges" << endl;

        ITstream& edgesStream(meshDescription.lookup("edges"));

        label nEdges = 0;

        token firstToken(edgesStream);

        if (firstToken.isLabel())
        {
            nEdges = firstToken.labelToken();
            edges_.setSize(nEdges);
        }
        else
        {
            edgesStream.putBack(firstToken);
        }

        // Read beginning of edges
        edgesStream.readBegin("edges");

        nEdges = 0;

        token lastToken(edgesStream);
        while
        (
            !(
                lastToken.isPunctuation()
                && lastToken.pToken() == token::END_LIST
            )
        )
        {
            if (edges_.size() <= nEdges)
            {
                edges_.setSize(nEdges + 1);
            }

            edgesStream.putBack(lastToken);

            edges_.set
            (
                nEdges,
                curvedEdge::New(tmpBlockPoints, edgesStream)
            );

            nEdges++;

            edgesStream >> lastToken;
        }
        edgesStream.putBack(lastToken);

        // Read end of edges
        edgesStream.readEnd("edges");
    }
    else
    {
        Info<< nl << "There are no non-linear edges" << endl;
    }


    Info<< nl << "Creating blocks" << endl;
    {
        ITstream& blockDescriptorStream(meshDescription.lookup("blocks"));

        // read number of blocks in mesh
        label nBlocks = 0;

        token firstToken(blockDescriptorStream);

        if (firstToken.isLabel())
        {
            nBlocks = firstToken.labelToken();
            blocks.setSize(nBlocks);
        }
        else
        {
            blockDescriptorStream.putBack(firstToken);
        }

        // Read beginning of blocks
        blockDescriptorStream.readBegin("blocks");

        nBlocks = 0;

        token lastToken(blockDescriptorStream);
        while
        (
            !(
                lastToken.isPunctuation()
                && lastToken.pToken() == token::END_LIST
            )
        )
        {
            if (blocks.size() <= nBlocks)
            {
                blocks.setSize(nBlocks + 1);
            }

            blockDescriptorStream.putBack(lastToken);

            blocks.set
            (
                nBlocks,
                new block
                (
                    blockDescriptor
                    (
                        tmpBlockPoints,
                        edges_,
                        blockDescriptorStream
                    )
                )
            );

            topologyOK = topologyOK && blockLabelsOK
            (
                nBlocks,
                tmpBlockPoints,
                blocks[nBlocks].blockDef().blockShape()
            );

            nBlocks++;

            blockDescriptorStream >> lastToken;
        }
        blockDescriptorStream.putBack(lastToken);

        // Read end of blocks
        blockDescriptorStream.readEnd("blocks");
    }


    polyMesh* blockMeshPtr = NULL;

    Info<< nl << "Creating patches" << endl;

    if (meshDescription.found("patches"))
    {
        Info<< nl << "Reading patches section" << endl;

        faceListList tmpBlocksPatches;
        wordList patchNames;
        wordList patchTypes;
        wordList nbrPatchNames;

        topologyOK = topologyOK && readPatches
        (
            meshDescription,
            tmpBlockPoints,
            tmpBlocksPatches,
            patchNames,
            patchTypes,
            nbrPatchNames
        );

        if (!topologyOK)
        {
            FatalErrorIn("blockMesh::createTopology(IOdictionary&)")
                << "Cannot create mesh due to errors in topology, exiting !"
                << nl << exit(FatalError);
        }


        Info<< nl << "Creating block mesh topology" << endl;

        PtrList<cellShape> tmpBlockCells(blocks.size());
        createCellShapes(tmpBlockPoints, tmpBlockCells);


        Info<< nl << "Reading physicalType from existing boundary file" << endl;

        wordList patchPhysicalTypes(tmpBlocksPatches.size());

        preservePatchTypes
        (
            meshDescription.time(),
            meshDescription.time().constant(),
            polyMesh::meshSubDir,
            patchNames,
            patchTypes,
            defaultPatchName,
            defaultPatchType,
            patchPhysicalTypes
        );

        blockMeshPtr = new polyMesh
        (
            IOobject
            (
                "blockMesh",
                meshDescription.time().constant(),
                meshDescription.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            tmpBlockPoints,
            tmpBlockCells,
            tmpBlocksPatches,
            patchNames,
            patchTypes,
            defaultPatchName,
            defaultPatchType,
            patchPhysicalTypes
        );
    }
    else if (meshDescription.found("boundary"))
    {
        faceListList tmpBlocksPatches;
        PtrList<dictionary> patchDicts;

        topologyOK = topologyOK && readBoundary
        (
            meshDescription,
            tmpBlockPoints,
            tmpBlocksPatches,
            patchDicts
        );

        if (!topologyOK)
        {
            FatalErrorIn("blockMesh::createTopology(IOdictionary&)")
                << "Cannot create mesh due to errors in topology, exiting !"
                << nl << exit(FatalError);
        }


        Info<< nl << "Creating block mesh topology" << endl;

        PtrList<cellShape> tmpBlockCells(blocks.size());
        createCellShapes(tmpBlockPoints, tmpBlockCells);


        blockMeshPtr = new polyMesh
        (
            IOobject
            (
                "blockMesh",
                meshDescription.time().constant(),
                meshDescription.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            tmpBlockPoints,
            tmpBlockCells,
            tmpBlocksPatches,
            patchDicts,
            defaultPatchName,
            defaultPatchType
        );
    }


    checkBlockMesh(*blockMeshPtr);

    return blockMeshPtr;
}


// ************************************************************************* //

