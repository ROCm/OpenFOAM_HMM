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
    Additional mesh accounting (Ensight)

\*---------------------------------------------------------------------------*/

PtrList<ensightCase> ensightCases(regionNames.size());
PtrList<ensightMesh> ensightMeshes(regionNames.size());

PtrList<faMesh> meshesFa(regionNames.size());
PtrList<ensightCase> ensightCasesFa(regionNames.size());
PtrList<ensightFaMesh> ensightMeshesFa(regionNames.size());

{
    forAll(regionNames, regioni)
    {
        const fvMesh& mesh = meshes[regioni];

        const word& regionName = regionNames[regioni];
        const word& regionDir =
        (
            regionName != polyMesh::defaultRegion
          ? regionName
          : word::null
        );

        fileName ensCasePath(outputDir);
        word ensCaseName(args.globalCaseName());

        if (!regionDir.empty())
        {
            ensCaseName = regionName;
            ensCasePath /= regionName;

            // Handle very rare naming collision with Ensight directories
            if (regionName == "data")
            {
                ensCasePath += "-region";
            }
        }

        ensightMeshes.set
        (
            regioni,
            new ensightMesh(mesh, writeOpts)
        );

        // New ensight case file, initialize header etc.
        ensightCases.set
        (
            regioni,
            new ensightCase(ensCasePath, ensCaseName, caseOpts)
        );

        if (doFiniteArea)
        {
            autoPtr<faMesh> faMeshPtr;

            const bool oldThrowingError = FatalError.throwing(true);
            try
            {
                faMeshPtr.reset(new faMesh(mesh));
            }
            catch (const Foam::error& err)
            {
                faMeshPtr.reset(nullptr);
            }
            FatalError.throwing(oldThrowingError);

            if (faMeshPtr)
            {
                ensightCasesFa.set
                (
                    regioni,
                    new ensightCase
                    (
                        ensCasePath/"finite-area",
                        "finite-area",
                        caseOpts
                    )
                );

                meshesFa.set(regioni, faMeshPtr);

                ensightMeshesFa.set
                (
                    regioni,
                    new ensightFaMesh(meshesFa[regioni])
                );
            }
        }
    }
}


// ************************************************************************* //
