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
    Info<< nl
        << "Writing edges in obj format" << endl;

    word outputName("finiteArea-edges.obj");

    if (Pstream::parRun())
    {
        outputName = word
        (
            "finiteArea-edges-proc"
          + Foam::name(Pstream::myProcNo())
          + ".obj"
        );
    }

    OBJstream os(runTime.globalPath()/outputName);

    os.writeQuoted
    (
        ("# " + outputName + "\n"),
        false
    );

    os.write(aMesh.patch().edges(), aMesh.patch().localPoints());
}


// ************************************************************************* //
