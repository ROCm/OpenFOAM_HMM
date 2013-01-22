scalar DiNum = -GREAT;

forAll(solidRegions, i)
{
    #include "setRegionSolidFields.H"

    tmp<volScalarField> magKappa;
    if (thermo.isotropic())
    {
        magKappa = thermo.kappa();
    }
    else
    {
        magKappa = mag(thermo.Kappa());
    }

    DiNum = max
    (
        solidRegionDiffNo
        (
            solidRegions[i],
            runTime,
            rho*cp,
            magKappa()
        ),
        DiNum
    );

}
