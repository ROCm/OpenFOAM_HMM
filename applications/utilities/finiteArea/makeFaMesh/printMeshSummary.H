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
    Summary of faMesh information

\*---------------------------------------------------------------------------*/

{
    const faBoundaryMesh& patches = aMesh.boundary();
    const label nNonProcessor = patches.nNonProcessor();
    const label nPatches = patches.size();

    Info<< "----------------" << nl
        << "Mesh Information" << nl
        << "----------------" << nl
        << "  " << "boundingBox: " << boundBox(aMesh.points()) << nl
        << "  " << "nFaces: " << returnReduce(aMesh.nFaces(), sumOp<label>())
        << nl;


    Info<< "----------------" << nl
        << "Patches" << nl
        << "----------------" << nl;

    for (label patchi = 0; patchi < nNonProcessor; ++patchi)
    {
        const faPatch& p = patches[patchi];

         Info<< "  " << "patch " << p.index()
             << " (size: " << returnReduce(p.size(), sumOp<label>())
             << ") name: " << p.name()
             << nl;
    }

    // Geometry information
    Info<< nl;
    {
        scalarMinMax limit(gMinMax(aMesh.S().field()));
        Info<< "Face area:" << nl
            << "    min = " << limit.min() << " max = "  << limit.max() << nl;
    }

    {
        scalarMinMax limit(minMax(aMesh.magLe().primitiveField()));

        // Include processor boundaries into 'internal' edges
        if (Pstream::parRun())
        {
            for (label patchi = nNonProcessor; patchi < nPatches; ++patchi)
            {
                limit.add(minMax(aMesh.magLe().boundaryField()[patchi]));
            }

            reduce(limit, minMaxOp<scalar>());
        }

        Info<< "Edge length (internal):" << nl
            << "    min = " << limit.min() << " max = "  << limit.max() << nl;


        // Include (non-processor) boundaries
        for (label patchi = 0; patchi < nNonProcessor; ++patchi)
        {
            limit.add(minMax(aMesh.magLe().boundaryField()[patchi]));
        }

        if (Pstream::parRun())
        {
            reduce(limit, minMaxOp<scalar>());
        }

        Info<< "Edge length:" << nl
            << "    min = " << limit.min()
            << " max = "  << limit.max() << nl;
    }

    // Not particularly meaningful
    #if 0
    {
        MinMax<vector> limit(gMinMax(aMesh.faceAreaNormals().field()));

        Info<< "Face area normals:" << nl
            << "    min = " << limit.min() << " max = " << limit.max() << nl;
    }
    #endif
}


// ************************************************************************* //
