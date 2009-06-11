{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "zero",
                dimensionSet(1, -3, -1, 0, 0, 0, 0),
                0.0
            )
        )
    );

    scalarField& dQ = tdQ();

    scalarField cp(dQ.size(), 0.0);

    forAll(Y, i)
    {
        volScalarField RRi = chemistry.RR(i);

        forAll(h, celli)
        {
            scalar Ti = T[celli];
            cp[celli] += Y[i][celli]*chemistry.specieThermo()[i].Cp(Ti);
            scalar hi = chemistry.specieThermo()[i].h(Ti);
            scalar RR = RRi[celli];
            dQ[celli] -= hi*RR;
        }
    }

    forAll(dQ, celli)
    {
        dQ[celli] /= cp[celli];
    }

    tdQ().write();
}
