{
    volScalarField& h = hSolid[i];

    h = allToProcMappers[procI].decomposeField(allh(), true);
    h.oldTime().timeIndex() = allh().oldTime().timeIndex();

    T += (h-h.oldTime())/cp;
    // Correct T boundary conditions and update h boundary
    // conditions accordingly.
    volScalarField::GeometricBoundaryField Told = T.boundaryField();
    T.correctBoundaryConditions();
    h.boundaryField() +=
        cp.boundaryField()
      * (T.boundaryField()-Told);
}
