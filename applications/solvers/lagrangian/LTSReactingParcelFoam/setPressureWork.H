DpDt == dimensionedScalar("zero", DpDt.dimensions(), 0.0);

if (pressureWork)
{
    surfaceScalarField phiU("phiU", phi/fvc::interpolate(rho));

    DpDt += fvc::div(phiU*fvc::interpolate(p)) - p*fvc::div(phiU);

    if (pressureWorkTimeDerivative)
    {
        DpDt += fvc::ddt(p);
    }
}
