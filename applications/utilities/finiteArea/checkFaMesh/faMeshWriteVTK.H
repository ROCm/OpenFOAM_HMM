/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    VTK output of faMesh with some geometric or debug fields

\*---------------------------------------------------------------------------*/

{
    // finiteArea - faces
    vtk::uindirectPatchWriter writer
    (
        aMesh.patch(),
        // vtk::formatType::INLINE_ASCII,
        fileName
        (
            aMesh.mesh().time().globalPath() / "finiteArea"
        )
    );

    writer.writeGeometry();

    globalIndex procAddr(aMesh.nFaces());
    labelList cellIDs;

    if (Pstream::master())
    {
        cellIDs.resize(procAddr.totalSize());
        for (const labelRange& range : procAddr.ranges())
        {
            auto slice = cellIDs.slice(range);
            slice = identity(range);
        }
    }

    // CellData
    writer.beginCellData(4);
    writer.writeProcIDs();
    writer.write("cellID", cellIDs);
    writer.write("area", aMesh.S().field());
    writer.write("normal", aMesh.faceAreaNormals());

    // PointData
    writer.beginPointData(1);
    writer.write("normal", aMesh.pointAreaNormals());

    Info<< nl
        << "Wrote faMesh in vtk format: " << writer.output().name() << nl;
}

{
    // finiteArea - edges
    vtk::lineWriter writer
    (
        aMesh.points(),
        aMesh.edges(),
        // vtk::formatType::INLINE_ASCII,
        fileName
        (
            aMesh.mesh().time().globalPath() / "finiteArea-edges"
        )
    );

    writer.writeGeometry();

    // CellData
    writer.beginCellData(4);
    writer.writeProcIDs();
    {
        // Use primitive patch order
        Field<scalar> fld
        (
            faMeshTools::flattenEdgeField(aMesh.magLe(), true)
        );
        writer.write("magLe", fld);
    }

    // PointData
    writer.beginPointData(1);
    writer.write("normal", aMesh.pointAreaNormals());

    Info<< nl
        << "Wrote faMesh in vtk format: " << writer.output().name() << nl;
}

{
    // finiteArea - edgeCentres
    // (no other convenient way to display vectors on the edges)
    vtk::lineWriter writer
    (
        aMesh.edgeCentres(),
        edgeList::null(),
        // vtk::formatType::INLINE_ASCII,
        fileName
        (
            aMesh.mesh().time().globalPath() / "finiteArea-edgesCentres"
        )
    );

    writer.writeGeometry();

    // PointData
    writer.beginPointData(4);
    {
        // Use primitive patch order
        Field<vector> fld
        (
            faMeshTools::flattenEdgeField(aMesh.Le(), true)
        );
        writer.write("Le", fld);
    }
    {
        // Use primitive patch order
        Field<vector> fld
        (
            faMeshTools::flattenEdgeField(aMesh.edgeAreaNormals(), true)
        );
        writer.write("normal", fld);
    }

    Info<< nl
        << "Wrote faMesh in vtk format: " << writer.output().name() << nl;
}


// ************************************************************************* //
