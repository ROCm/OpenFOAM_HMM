/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    VTK output of blockMesh topology blocks

\*---------------------------------------------------------------------------*/

{
    // Non-legacy and ASCII (mesh is small, want readable output)
    const vtk::outputOptions writeOpts = vtk::formatType::INLINE_ASCII;

    const polyMesh& topoMesh = blocks.topology();
    const vtk::vtuCells topoCells(topoMesh, writeOpts);

    vtk::internalMeshWriter writer
    (
        topoMesh,
        topoCells,
        writeOpts,
        runTime.path()/"blockTopology"
    );

    Info<< "Writing block topology in vtk format: "
        << args.relativePath(writer.output()).c_str() << endl;

    writer.writeGeometry();
    writer.beginCellData();
    writer.writeCellIDs();

    // No cell decomposition, so there is a 1-to-1 correspondence between
    // input blocks and VTK output cells.

    vectorField localNormal(blocks.size());

    // Generate local normals as fields for visualisation
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        const label faceMin = label(2*cmpt);
        const label faceMax = faceMin+1;

        localNormal.resize(blocks.size());

        forAll(blocks, blocki)
        {
            const cellShape& shape = blocks[blocki].blockShape();
            const pointField& verts = blocks[blocki].vertices();

            if (shape.model() == cellModel::ref(cellModel::HEX))
            {
                localNormal[blocki] =
                (
                    // 30% cell-width as arbitrary value for vector length
                    0.3*mag
                    (
                        shape.face(faceMax).centre(verts)
                      - shape.face(faceMin).centre(verts)
                    )
                  * normalised
                    (
                        // Weigh by area to avoid influence of zero-faces
                        shape.face(faceMax).areaNormal(verts)
                      - shape.face(faceMin).areaNormal(verts)
                    )
                );
            }
            else
            {
                // Could handle other shapes (how?)
                localNormal[blocki] = Zero;
            }
        }

        // Additional safety (should not happen)
        localNormal.resize(topoMesh.nCells(), Zero);

        writer.writeCellData
        (
            word("local-direction" + name(cmpt)),
            localNormal
        );
    }

    // if (topoMesh.nCells() != blocks.size())
    // {
    //     Info<< "Warning: indicated normals may be incorrect" << nl;
    // }
}


// ************************************************************************* //
