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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::polyMesh* Foam::blockMesh::createTopology(IOdictionary& dict)
{
    bool topologyOK = true;

    blockMesh& blocks = *this;

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
        Info<< "Creating curved edges" << endl;

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
    else
    {
        Info<< "No non-linear edges defined" << endl;
    }


    //
    // Create the blocks
    //
    Info<< "Creating topology blocks" << endl;
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


    //
    // Create the patches
    //
    Info<< "Creating topology patches" << endl;

    faceListList tmpBlocksPatches;
    wordList patchNames;
    wordList patchTypes;

    {
        ITstream& is(dict.lookup("patches"));

        // read number of patches in mesh
        label nPatches = 0;

        token firstToken(is);

        if (firstToken.isLabel())
        {
            nPatches = firstToken.labelToken();

            tmpBlocksPatches.setSize(nPatches);
            patchNames.setSize(nPatches);
            patchTypes.setSize(nPatches);
        }
        else
        {
            is.putBack(firstToken);
        }

        // Read beginning of blocks
        is.readBegin("patches");

        nPatches = 0;

        token lastToken(is);
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
            }

            is.putBack(lastToken);

            is
                >> patchTypes[nPatches]
                >> patchNames[nPatches]
                >> tmpBlocksPatches[nPatches];


            // Catch multiple patches asap.
            for (label i = 0; i < nPatches; i++)
            {
                if (patchNames[nPatches] == patchNames[i])
                {
                    FatalErrorIn
                    (
                        "blockMesh::createTopology(IOdictionary&)"
                    )   << "Duplicate patch " << patchNames[nPatches]
                        << " at line " << is.lineNumber()
                        << ". Exiting !" << nl
                        << exit(FatalError);
                }
            }

            topologyOK = topologyOK && patchLabelsOK
            (
                nPatches,
                blockPointField_,
                tmpBlocksPatches[nPatches]
            );

            nPatches++;

            is >> lastToken;
        }
        is.putBack(lastToken);

        // Read end of blocks
        is.readEnd("patches");
    }


    if (!topologyOK)
    {
        FatalErrorIn("blockMesh::createTopology(IOdictionary&)")
            << "Cannot create mesh due to errors in topology, exiting !" << nl
            << exit(FatalError);
    }


    //
    // Create the topology
    //
    Info<< "Creating topology mesh" << endl;

    PtrList<cellShape> tmpBlockShapes(blocks.size());
    forAll(blocks, blockI)
    {
        tmpBlockShapes.set
        (
            blockI,
            new cellShape(blocks[blockI].blockShape())
        );

        if (tmpBlockShapes[blockI].mag(blockPointField_) < 0.0)
        {
            WarningIn
            (
                "blockMesh::createTopology(IOdictionary&)"
            )   << "negative volume block : " << blockI
                << ", probably defined inside-out" << endl;
        }
    }

    wordList patchPhysicalTypes(tmpBlocksPatches.size());

    preservePatchTypes
    (
        dict.time(),
        dict.time().constant(),
        polyMesh::meshSubDir,
        patchNames,
        patchTypes,
        defaultPatchName,
        defaultPatchType,
        patchPhysicalTypes
    );


    // construct the topology as its own mesh
    polyMesh* blockMeshPtr = new polyMesh
    (
        IOobject
        (
            "blockMesh",
            dict.time().constant(),
            dict.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        xferCopy(blockPointField_),   // copy these points, do NOT move
        tmpBlockShapes,
        tmpBlocksPatches,
        patchNames,
        patchTypes,
        defaultPatchName,
        defaultPatchType,
        patchPhysicalTypes
    );

    checkBlockMesh(*blockMeshPtr);

    return blockMeshPtr;
}


// ************************************************************************* //
