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
    Summary of faMesh information

\*---------------------------------------------------------------------------*/

{
    const faBoundaryMesh& patches = aMesh.boundary();
    const label nNonProcessor = patches.nNonProcessor();
    const label nPatches = patches.size();

    label nLocalProcEdges = 0;
    if (Pstream::parRun())
    {
        for (const faPatch& fap : patches)
        {
            const auto* cpp = isA<processorFaPatch>(fap);

            if (cpp)
            {
                nLocalProcEdges += fap.nEdges();
            }
        }
    }

    const labelList nFaces
    (
        UPstream::listGatherValues<label>(aMesh.nFaces())
    );
    const labelList nPoints
    (
        UPstream::listGatherValues<label>(aMesh.nPoints())
    );
    const labelList nEdges
    (
        UPstream::listGatherValues<label>(aMesh.nEdges())
    );
    const labelList nIntEdges
    (
        UPstream::listGatherValues<label>(aMesh.nInternalEdges())
    );

    // The "real" (non-processor) boundary edges
    const labelList nBndEdges
    (
        UPstream::listGatherValues<label>
        (
            aMesh.nBoundaryEdges() - nLocalProcEdges
        )
    );
    const labelList nProcEdges
    (
        UPstream::listGatherValues<label>(nLocalProcEdges)
    );


    // Format output as
    //  Number of faces: ...
    //         per-proc: (...)

    const auto reporter =
        [&](const char* tag, const labelList& list)
        {
            Info<< "  Number of " << tag << ": " << sum(list) << nl;
            if (Pstream::parRun())
            {
                int padding = static_cast<int>
                (
                    // strlen("  Number of ") - strlen("per-proc")
                    (12 - 8)
                  + strlen(tag)
                );

                do { Info<< ' '; } while (--padding > 0);

                Info<< "per-proc: " << flatOutput(list) << nl;
            }
        };


    Info<< "----------------" << nl
        << "Mesh Information" << nl
        << "----------------" << nl
        << "  " << "boundingBox: " << boundBox(aMesh.points()) << nl;

    if (Pstream::master())
    {
        reporter("faces", nFaces);
        reporter("points", nPoints);
        reporter("edges", nEdges);
        reporter("internal edges", nIntEdges);
        reporter("boundary edges", nBndEdges);

        if (Pstream::parRun())
        {
            reporter("processor edges", nProcEdges);
        }
    }

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
            << "    min = " << limit.min() << " max = "  << limit.max() << nl;
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
