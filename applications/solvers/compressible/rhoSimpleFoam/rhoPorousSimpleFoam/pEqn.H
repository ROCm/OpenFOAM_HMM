{
    rho = thermo.rho();
    rho = max(rho, rhoMin);
    rho = min(rho, rhoMax);
    rho.relax();

    const volScalarField& psi = thermo.psi();

    volVectorField HbyA("HbyA", U);

    if (pressureImplicitPorosity)
    {
        HbyA = trTU() & (UEqn() == sources(rho, U))().H();
    }
    else
    {
        HbyA = trAU()*(UEqn() == sources(rho, U))().H();
    }

    UEqn.clear();

    bool closedVolume = false;

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho*HbyA) & mesh.Sf()
    );

    mrfZones.relativeFlux(fvc::interpolate(rho), phiHbyA);

    closedVolume = adjustPhi(phiHbyA, U, p);

    while (simple.correctNonOrthogonal())
    {
        tmp<fvScalarMatrix> tpEqn;

        if (pressureImplicitPorosity)
        {
            tpEqn =
            (
                fvm::laplacian(rho*trTU(), p)
              + sources(psi, p, rho.name())
             ==
                fvc::div(phiHbyA)
            );
        }
        else
        {
            tpEqn =
            (
                fvm::laplacian(rho*trAU(), p)
              + sources(psi, p, rho.name())
             ==
                fvc::div(phiHbyA)
            );
        }

        tpEqn().setReference(pRefCell, pRefValue);

        sources.constrain(tpEqn(), rho.name());

        tpEqn().solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - tpEqn().flux();
        }
    }

    #include "incompressible/continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    if (pressureImplicitPorosity)
    {
        U = HbyA - (trTU() & fvc::grad(p));
    }
    else
    {
        U = HbyA - trAU()*fvc::grad(p);
    }

    U.correctBoundaryConditions();
    sources.correct(U);

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
    }

    rho = thermo.rho();
    rho = max(rho, rhoMin);
    rho = min(rho, rhoMax);
    rho.relax();
    Info<< "rho max/min : "
        << max(rho).value() << " "
        << min(rho).value() << endl;
}
