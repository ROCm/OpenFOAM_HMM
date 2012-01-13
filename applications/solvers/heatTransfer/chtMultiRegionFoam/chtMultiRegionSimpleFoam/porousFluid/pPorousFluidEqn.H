porousRho = porousThermo.rho();
porousRho = max(porousRho, rhoMin);
porousRho = min(porousRho, rhoMax);
porousRho.relax();

volScalarField rAUPorous(1.0/porousUEqn().A());

porousU = rAUPorous*porousUEqn().H();
porousUEqn.clear();

bool closedVolume = false;

porousPhi =
    fvc::interpolate(porousRho)
    *(fvc::interpolate(porousU) & porousMesh.Sf());

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{

    fvScalarMatrix pEqn
    (
        fvm::laplacian(porousRho*rAUPorous, porousP) ==  fvc::div(porousPhi)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (nonOrth == nNonOrthCorr)
    {
        porousPhi -= pEqn.flux();

    }
}

porousP.relax();

porousU -= rAUPorous*fvc::grad(porousP);
porousU.correctBoundaryConditions();

if (closedVolume)
{
    porousP += (initialMass - fvc::domainIntegrate(porousPsi*porousP))
        /fvc::domainIntegrate(porousPsi);
}

porousRho = porousThermo.rho();
porousRho = max(porousRho, rhoMin);
porousRho = min(porousRho, rhoMax);
porousRho.relax();

Info<< "rho max/min : "
    << max(porousRho).value() << " "
    << min(porousRho).value() << endl;
