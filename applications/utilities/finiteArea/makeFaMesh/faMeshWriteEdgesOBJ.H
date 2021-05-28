/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    OBJ output of faMesh edges

\*---------------------------------------------------------------------------*/

{
    Info<< "Writing edges in obj format" << endl;

    word outputName("faMesh-edges.obj");

    if (Pstream::parRun())
    {
        outputName = word
        (
            "faMesh-edges-" + Foam::name(Pstream::myProcNo()) + ".obj"
        );
    }

    OBJstream os(runTime.globalPath()/outputName);

    os.writeQuoted
    (
        ("# " + outputName + "\n"),
        false
    );

    os.write(areaMesh.patch().edges(), areaMesh.patch().localPoints());
}


// ************************************************************************* //
