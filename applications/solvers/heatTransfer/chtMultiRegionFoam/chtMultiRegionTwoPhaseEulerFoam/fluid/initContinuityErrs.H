PtrList<uniformDimensionedScalarField> cumulativeContErrIO(fluidRegions.size());
forAll(cumulativeContErrIO, i)
{
    #include "setRegionFluidFields.H"
    cumulativeContErrIO.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "cumulativeContErr",
                runTime.timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            dimensionedScalar(dimless, Zero)
        )
    );
}

UPtrList<scalar> cumulativeContErr(cumulativeContErrIO.size());
forAll(cumulativeContErrIO, i)
{
    cumulativeContErr.set(i, &cumulativeContErrIO[i].value());
}
