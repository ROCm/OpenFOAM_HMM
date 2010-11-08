DpDt == dimensionedScalar("zero", DpDt.dimensions(), 0.0);

if (eWork)
{
    DpDt += -p*fvc::div(phi/fvc::interpolate(rho));
}
if (hWork)
{
    DpDt += fvc::div(phi/fvc::interpolate(rho)*fvc::interpolate(p));
}
