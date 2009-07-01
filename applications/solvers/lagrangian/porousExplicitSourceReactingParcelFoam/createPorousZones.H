    Info<< "Creating porous zones" << nl << endl;

    porousZones pZones(mesh);
    Switch pressureImplicitPorosity(false);

    label nUCorr = 0;
    if (pZones.size())
    {
        // nUCorrectors for pressureImplicitPorosity
        if (mesh.solutionDict().subDict("PISO").found("nUCorrectors"))
        {
            mesh.solutionDict().subDict("PISO").lookup("nUCorrectors")
                >> nUCorr;
        }

        if (nUCorr > 0)
        {
            pressureImplicitPorosity = true;
        }
    }

