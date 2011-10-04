{
    volScalarField& h = hSolid[i];

    procPhi.setSize(nAllRegions);
    procAlpha.setSize(nAllRegions);

    rmap
    (
        allRho(),
        rho,
        faceProcAddressing[procI],
        cellProcAddressing[procI],
        boundaryProcAddressing[procI]
    );

    // Necessary? Probably only for boundary values since bcs on
    // h are not the same as those on allh
    rmap
    (
        allh(),
        h,
        faceProcAddressing[procI],
        cellProcAddressing[procI],
        boundaryProcAddressing[procI]
    );


    tmp<volScalarField> Kcp(K/cp);

    rmap
    (
        allVolAlpha(),
        Kcp(),
        faceProcAddressing[procI],
        cellProcAddressing[procI],
        boundaryProcAddressing[procI]
    );

    procAlpha.set(procI, fvc::interpolate(Kcp));

    // allSource is initialised to zero already
    //rmap
    //(
    //    allSource(),
    //    volScalarField
    //    (
    //        IOobject
    //        (
    //            "procSource",
    //            runTime.timeName(),
    //            mesh,
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE
    //        ),
    //        mesh,
    //        dimensionedScalar
    //        (
    //            "procSource",
    //            allh().dimensions()*dimDensity/dimTime,
    //            0.0
    //        )
    //    ),
    //    faceProcAddressing[procI],
    //    cellProcAddressing[procI],
    //    boundaryProcAddressing[procI]
    //);

    procPhi.set
    (
        procI,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "phi",
                dimDensity*dimVelocity*dimArea,
                0.0
            )
        )
    );
}
