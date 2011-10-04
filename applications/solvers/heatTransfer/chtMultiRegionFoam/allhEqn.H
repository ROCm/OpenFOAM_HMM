
// Get mapped alpha (surfaceScalarField)
tmp<surfaceScalarField> tallAlpha = procToAllMapper().reconstructFvSurfaceField
    (
        IOobject
        (
            "alpha",
            allMesh().time().timeName(),
            allMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        procAlpha
    );

// Get alpha from harmonic interpolation of vol quantities
// (Note: really only needed at internal faces originating patches
//  inbetween regions)
tmp<surfaceScalarField> allHarmonicAlpha
(
    harmonic(allMesh()).interpolate(allVolAlpha())
);

// Loop over all fluid and solid regions to transfer
// allHarmonicAlpha to allAlpha
surfaceScalarField& allAlpha = tallAlpha();
forAll(boundaryProcAddressing, procI)
{
    forAll(boundaryProcAddressing[procI], patchI)
    {
        if (boundaryProcAddressing[procI][patchI] == -1)
        {
            // Interface patch
            const labelList::subList cp =
                procMeshes[procI].boundary()[patchI].patchSlice
                (
                    faceProcAddressing[procI]
                );

            forAll(cp, faceI)
            {
                label curF = mag(cp[faceI])-1;
                if (curF < allMesh().nInternalFaces())
                {
                    allAlpha[curF] = allHarmonicAlpha()[curF];
                }
            }
        }
    }
}


tmp<surfaceScalarField> allPhi
(
    procToAllMapper().reconstructFvSurfaceField
    (
        IOobject
        (
            "phi",
            allMesh().time().timeName(),
            allMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        procPhi
    )
);

// So we have nNonOrthCorr
//#include "readSolidMultiRegionPIMPLEControls.H"
//for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix hEqn
    (
        fvm::ddt(allRho(), allh())
      + fvm::div(allPhi(), allh())
      - fvm::laplacian(allAlpha, allh())
     ==
        allSource()
    );

    hEqn.relax();
    hEqn.solve(allMesh().solver(allh().select(finalIter)));
}
