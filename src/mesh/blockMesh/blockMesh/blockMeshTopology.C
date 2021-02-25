/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "blockMesh.H"
#include "blockMeshTools.H"
#include "Time.H"
#include "preservePatchTypes.H"
#include "emptyPolyPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Source>
void Foam::blockMesh::checkPatchLabels
(
    const Source& source,
    const word& patchName,
    const pointField& points,
    faceList& patchFaces
) const
{
    forAll(patchFaces, facei)
    {
        face& f = patchFaces[facei];

        // Replace (<block> <face>) face description
        // with the corresponding block face
        if (f.size() == 2)
        {
            const label bi = f[0];
            const label fi = f[1];

            if (bi >= size())
            {
                FatalIOErrorInFunction(source)
                    << "Block index out of range for patch face " << f << nl
                    << "    Number of blocks = " << size()
                    << ", index = " << f[0] << nl
                    << "    on patch " << patchName << ", face " << facei
                    << exit(FatalIOError);
            }
            else if (fi >= operator[](bi).blockShape().faces().size())
            {
                FatalIOErrorInFunction(source)
                    << "Block face index out of range for patch face " << f
                    << nl
                    << "    Number of block faces = "
                    << operator[](bi).blockShape().faces().size()
                    << ", index = " << f[1] << nl
                    << "    on patch " << patchName << ", face " << facei
                    << exit(FatalIOError);
            }
            else
            {
                f = operator[](bi).blockShape().faces()[fi];
            }
        }
        else
        {
            for (const label pointi : f)
            {
                if (pointi < 0 || pointi >= points.size())
                {
                    FatalIOErrorInFunction(source)
                        << "Point label " << pointi
                        << " out of range 0.." << points.size() - 1 << nl
                        << "    on patch " << patchName << ", face " << facei
                        << exit(FatalIOError);
                }
            }
        }
    }
}


void Foam::blockMesh::readPatches
(
    const dictionary& meshDescription,
    faceListList& tmpBlocksPatches,
    wordList& patchNames,
    wordList& patchTypes,
    wordList& nbrPatchNames
)
{
    // Collect all variables
    dictionary varDict(meshDescription.subOrEmptyDict("namedVertices"));
    varDict.merge(meshDescription.subOrEmptyDict("namedBlocks"));


    ITstream& patchStream(meshDescription.lookup("patches"));

    // Read number of patches in mesh
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
    while (!lastToken.isPunctuation(token::END_LIST))
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

        // Read patch faces
        tmpBlocksPatches[nPatches] = blockMeshTools::read<face>
        (
            patchStream,
            varDict
        );


        // Check for multiple patches
        for (label i = 0; i < nPatches; i++)
        {
            if (patchNames[nPatches] == patchNames[i])
            {
                FatalIOErrorInFunction(patchStream)
                    << "Duplicate patch " << patchNames[nPatches]
                    << " at line " << patchStream.lineNumber()
                    << exit(FatalIOError);
            }
        }

        checkPatchLabels
        (
            patchStream,
            patchNames[nPatches],
            vertices_,
            tmpBlocksPatches[nPatches]
        );

        nPatches++;


        // Split old style cyclics

        if (patchTypes[nPatches-1] == cyclicPolyPatch::typeName)
        {
            word halfA = patchNames[nPatches-1] + "_half0";
            word halfB = patchNames[nPatches-1] + "_half1";

            FatalIOErrorInFunction(patchStream)
                << "Old-style cyclic definition."
                << " Splitting patch "
                << patchNames[nPatches-1] << " into two halves "
                << halfA << " and " << halfB << endl
                << "    Alternatively use new 'boundary' dictionary syntax."
                << exit(FatalIOError);

            // Add extra patch
            if (tmpBlocksPatches.size() <= nPatches)
            {
                tmpBlocksPatches.setSize(nPatches + 1);
                patchNames.setSize(nPatches + 1);
                patchTypes.setSize(nPatches + 1);
                nbrPatchNames.setSize(nPatches + 1);
            }

            // Update halfA info
            patchNames[nPatches-1] = halfA;
            nbrPatchNames[nPatches-1] = halfB;

            // Update halfB info
            patchTypes[nPatches] = patchTypes[nPatches-1];
            patchNames[nPatches] = halfB;
            nbrPatchNames[nPatches] = halfA;

            // Split faces
            if ((tmpBlocksPatches[nPatches-1].size() % 2) != 0)
            {
                FatalIOErrorInFunction(patchStream)
                    << "Size of cyclic faces is not a multiple of 2 :"
                    << tmpBlocksPatches[nPatches-1]
                    << exit(FatalIOError);
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

        patchStream >> lastToken;
    }
    patchStream.putBack(lastToken);

    // Read end of blocks
    patchStream.readEnd("patches");
}


void Foam::blockMesh::readBoundary
(
    const dictionary& meshDescription,
    wordList& patchNames,
    faceListList& tmpBlocksPatches,
    PtrList<dictionary>& patchDicts
)
{
    // Collect all variables
    dictionary varDict(meshDescription.subOrEmptyDict("namedVertices"));
    varDict.merge(meshDescription.subOrEmptyDict("namedBlocks"));


    // Read like boundary file
    const PtrList<entry> patchesInfo
    (
        meshDescription.lookup("boundary")
    );

    patchNames.setSize(patchesInfo.size());
    tmpBlocksPatches.setSize(patchesInfo.size());
    patchDicts.setSize(patchesInfo.size());

    forAll(tmpBlocksPatches, patchi)
    {
        const entry& patchInfo = patchesInfo[patchi];

        if (!patchInfo.isDict())
        {
            FatalIOErrorInFunction(meshDescription)
                << "Entry " << patchInfo << " in boundary section is not a"
                << " valid dictionary."
                << exit(FatalIOError);
        }

        patchNames[patchi] = patchInfo.keyword();

        // Construct patch dictionary
        patchDicts.set(patchi, new dictionary(patchInfo.dict()));

        // Read block faces
        tmpBlocksPatches[patchi] = blockMeshTools::read<face>
        (
            patchDicts[patchi].lookup("faces"),
            varDict
        );


        checkPatchLabels
        (
            patchInfo.dict(),
            patchNames[patchi],
            vertices_,
            tmpBlocksPatches[patchi]
        );
    }
}


void Foam::blockMesh::createCellShapes
(
    cellShapeList& tmpBlockCells
)
{
    const blockMesh& blocks = *this;

    tmpBlockCells.setSize(blocks.size());
    forAll(blocks, blocki)
    {
        tmpBlockCells[blocki] = blocks[blocki].blockShape();
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::polyMesh> Foam::blockMesh::createTopology
(
    const IOdictionary& meshDescription,
    const word& regionName
)
{
    blockList& blocks = *this;

    word defaultPatchName = "defaultFaces";
    word defaultPatchType = emptyPolyPatch::typeName;

    // Read the names/types for the unassigned patch faces
    // this is a bit heavy handed (and ugly), but there is currently
    // no easy way to rename polyMesh patches subsequently
    if (const dictionary* dictPtr = meshDescription.findDict("defaultPatch"))
    {
        dictPtr->readIfPresent("name", defaultPatchName);
        dictPtr->readIfPresent("type", defaultPatchType);
    }

    // Optional 'scale' factor. Was 'convertToMeters' until OCT-2008
    meshDescription.readIfPresentCompat
    (
        "scale",
        {{"convertToMeters", 1012}},  // Mark as changed from 2010 onwards
        scaleFactor_
    );

    // Treat (scale <= 0) as scaling == 1 (no scaling).
    if (scaleFactor_ <= 0)
    {
        scaleFactor_ = 1.0;
    }

    // Read the block edges
    if (meshDescription.found("edges"))
    {
        if (verbose_)
        {
            Info<< "Creating block edges" << endl;
        }

        blockEdgeList edges
        (
            meshDescription.lookup("edges"),
            blockEdge::iNew(meshDescription, geometry_, vertices_)
        );

        edges_.transfer(edges);
    }
    else if (verbose_)
    {
        Info<< "No non-linear block edges defined" << endl;
    }


    // Read the block faces
    if (meshDescription.found("faces"))
    {
        if (verbose_)
        {
            Info<< "Creating block faces" << endl;
        }

        blockFaceList faces
        (
            meshDescription.lookup("faces"),
            blockFace::iNew(meshDescription, geometry_)
        );

        faces_.transfer(faces);
    }
    else if (verbose_)
    {
        Info<< "No non-planar block faces defined" << endl;
    }


    // Create the blocks
    if (verbose_)
    {
        Info<< "Creating topology blocks" << endl;
    }
    {
        blockList blocks
        (
            meshDescription.lookup("blocks"),
            block::iNew(meshDescription, vertices_, edges_, faces_)
        );

        transfer(blocks);
    }


    autoPtr<polyMesh> blockMeshPtr;

    // Create the patches

    if (verbose_)
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

        readPatches
        (
            meshDescription,
            tmpBlocksPatches,
            patchNames,
            patchTypes,
            nbrPatchNames
        );

        Info<< nl << "Creating block mesh topology" << endl;

        cellShapeList tmpBlockCells(blocks.size());
        createCellShapes(tmpBlockCells);


        Info<< nl << "Reading physicalType from existing boundary file" << endl;

        PtrList<dictionary> patchDicts(patchNames.size());
        word defaultFacesType;

        preservePatchTypes
        (
            meshDescription.time(),
            meshDescription.time().constant(),
            polyMesh::meshSubDir,
            patchNames,
            patchDicts,
            defaultPatchName,
            defaultPatchType
        );


        // Add cyclic info (might not be present from older file)
        forAll(patchDicts, patchi)
        {
            if (!patchDicts.set(patchi))
            {
                patchDicts.set(patchi, new dictionary());
            }

            dictionary& dict = patchDicts[patchi];

            // Add but not override type
            if (!dict.found("type"))
            {
                dict.add("type", patchTypes[patchi], false);
            }
            else if (dict.get<word>("type") != patchTypes[patchi])
            {
                FatalIOErrorInFunction(meshDescription)
                    << "For patch " << patchNames[patchi]
                    << " overriding type '" << patchTypes[patchi]
                    << "' with '" << dict.get<word>("type")
                    << "' (read from boundary file)"
                    << exit(FatalIOError);
            }

            // Override neighbourpatch name
            if (nbrPatchNames[patchi] != word::null)
            {
                dict.set("neighbourPatch", nbrPatchNames[patchi]);
            }
        }

        blockMeshPtr = autoPtr<polyMesh>::New
        (
            IOobject
            (
                regionName,
                meshDescription.time().constant(),
                meshDescription.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointField(vertices_),   // Copy these points, do NOT move
            tmpBlockCells,
            tmpBlocksPatches,
            patchNames,
            patchDicts,
            defaultPatchName,
            defaultPatchType
        );
    }
    else if (meshDescription.found("boundary"))
    {
        wordList patchNames;
        faceListList tmpBlocksPatches;
        PtrList<dictionary> patchDicts;

        readBoundary
        (
            meshDescription,
            patchNames,
            tmpBlocksPatches,
            patchDicts
        );

        Info<< nl << "Creating block mesh topology" << endl;

        cellShapeList tmpBlockCells(blocks.size());
        createCellShapes(tmpBlockCells);

        blockMeshPtr = autoPtr<polyMesh>::New
        (
            IOobject
            (
                regionName,
                meshDescription.time().constant(),
                meshDescription.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointField(vertices_),   // Copy these points, do NOT move
            tmpBlockCells,
            tmpBlocksPatches,
            patchNames,
            patchDicts,
            defaultPatchName,
            defaultPatchType
        );
    }

    check(*blockMeshPtr, meshDescription);

    return blockMeshPtr;
}


// ************************************************************************* //
