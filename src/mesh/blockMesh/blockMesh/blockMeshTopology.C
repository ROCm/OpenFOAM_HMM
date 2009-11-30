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

#include "blockMesh.H"
#include "Time.H"
#include "preservePatchTypes.H"
#include "emptyPolyPatch.H"
#include "cyclicPolyPatch.H"


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

Foam::polyMesh* Foam::blockMesh::createTopology(IOdictionary& dict)
{
    bool topologyOK = true;

    blockList& blocks = *this;

    word defaultPatchName = "defaultFaces";
    word defaultPatchType = emptyPolyPatch::typeName;

    // get names/types for the unassigned patch faces
    // this is a bit heavy handed (and ugly), but there is currently
    // no easy way to rename polyMesh patches subsequently
    if (const dictionary* dictPtr = dict.subDictPtr("defaultPatch"))
    {
        dictPtr->readIfPresent("name", defaultPatchName);
        dictPtr->readIfPresent("type", defaultPatchType);
    }

    // optional 'convertToMeters' or 'scale'  scaling factor
    if (!dict.readIfPresent("convertToMeters", scaleFactor_))
    {
        dict.readIfPresent("scale", scaleFactor_);
    }


    //
    // get the non-linear edges in mesh
    //
    if (dict.found("edges"))
    {
        if (verboseOutput)
        {
            Info<< "Creating curved edges" << endl;
        }

        ITstream& is(dict.lookup("edges"));

        // read number of edges in mesh
        label nEdges = 0;

        token firstToken(is);

        if (firstToken.isLabel())
        {
            nEdges = firstToken.labelToken();
            edges_.setSize(nEdges);
        }
        else
        {
            is.putBack(firstToken);
        }

        // Read beginning of edges
        is.readBegin("edges");

        nEdges = 0;

        token lastToken(is);
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

            is.putBack(lastToken);

            edges_.set
            (
                nEdges,
                curvedEdge::New(blockPointField_, is)
            );

            nEdges++;

            is >> lastToken;
        }
        is.putBack(lastToken);

        // Read end of edges
        is.readEnd("edges");
    }
    else if (verboseOutput)
    {
        Info<< "No non-linear edges defined" << endl;
    }


    //
    // Create the blocks
    //
    if (verboseOutput)
    {
        Info<< "Creating topology blocks" << endl;
    }

    {
        ITstream& is(dict.lookup("blocks"));

        // read number of blocks in mesh
        label nBlocks = 0;

        token firstToken(is);

        if (firstToken.isLabel())
        {
            nBlocks = firstToken.labelToken();
            blocks.setSize(nBlocks);
        }
        else
        {
            is.putBack(firstToken);
        }

        // Read beginning of blocks
        is.readBegin("blocks");

        nBlocks = 0;

        token lastToken(is);
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

            is.putBack(lastToken);

            blocks.set
            (
                nBlocks,
                new block
                (
                    blockPointField_,
                    edges_,
                    is
                )
            );

            topologyOK = topologyOK && blockLabelsOK
            (
                nBlocks,
                blockPointField_,
                blocks[nBlocks].blockShape()
            );

            nBlocks++;

            is >> lastToken;
        }
        is.putBack(lastToken);

        // Read end of blocks
        is.readEnd("blocks");
    }


    polyMesh* blockMeshPtr = NULL;

    //
    // Create the patches
    //
    if (verboseOutput)
    {
        Info<< "Creating topology patches" << endl;
    }

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
            xferMove(tmpBlockPoints),
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
            xferMove(tmpBlockPoints),
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
