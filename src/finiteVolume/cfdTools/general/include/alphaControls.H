const dictionary& alphaControls = mesh.solverDict(alpha1.name());
label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));
label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
Switch MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));

if (nAlphaSubCycles > 1 && pimple.nCorrPIMPLE() != 1)
{
    FatalErrorIn(args.executable())
        << "Sub-cycling alpha is only allowed for PISO operation, "
           "i.e. when the number of outer-correctors = 1"
        << exit(FatalError);
}
