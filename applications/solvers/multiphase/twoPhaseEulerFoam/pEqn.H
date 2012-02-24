{
    surfaceScalarField alphaf(fvc::interpolate(alpha));
    surfaceScalarField betaf(scalar(1) - alphaf);

    volScalarField rUaA(1.0/UaEqn.A());
    volScalarField rUbA(1.0/UbEqn.A());

    phia.boundaryField() ==
        (fvc::interpolate(Ua) & mesh.Sf())().boundaryField();
    phib.boundaryField() ==
        (fvc::interpolate(Ub) & mesh.Sf())().boundaryField();

    rUaAf = fvc::interpolate(rUaA);
    surfaceScalarField rUbAf(fvc::interpolate(rUbA));

    volVectorField HbyAa("HbyAa", Ua);
    HbyAa = rUaA*UaEqn.H();

    volVectorField HbyAb("HbyAb", Ub);
    HbyAb = rUbA*UbEqn.H();

    surfaceScalarField phiDraga
    (
        fvc::interpolate(beta/rhoa*K*rUaA)*phib + rUaAf*(g & mesh.Sf())
    );

    if (g0.value() > 0.0)
    {
        phiDraga -= ppMagf*fvc::snGrad(alpha)*mesh.magSf();
    }

    if (kineticTheory.on())
    {
        phiDraga -= rUaAf*fvc::snGrad(kineticTheory.pa()/rhoa)*mesh.magSf();
    }

    surfaceScalarField phiDragb
    (
        fvc::interpolate(alpha/rhob*K*rUbA)*phia + rUbAf*(g & mesh.Sf())
    );

    // Fix for gravity on outlet boundary.
    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
    }

    surfaceScalarField phiHbyAa
    (
        (fvc::interpolate(HbyAa) & mesh.Sf())
      + fvc::ddtPhiCorr(rUaA, Ua, phia)
      + phiDraga
    );

    surfaceScalarField phiHbyAb
    (
        (fvc::interpolate(HbyAb) & mesh.Sf())
      + fvc::ddtPhiCorr(rUbA, Ub, phib)
      + phiDragb
    );

    surfaceScalarField phiHbyA("phiHbyA", alphaf*phiHbyAa + betaf*phiHbyAb);

    surfaceScalarField Dp
    (
        "Dp",
        alphaf*rUaAf/rhoa + betaf*rUbAf/rhob
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField SfGradp(pEqn.flux()/Dp);

            phia = phiHbyAa - rUaAf*SfGradp/rhoa;
            phib = phiHbyAb - rUbAf*SfGradp/rhob;
            phi = alphaf*phia + betaf*phib;

            p.relax();
            SfGradp = pEqn.flux()/Dp;

            Ua = HbyAa + fvc::reconstruct(phiDraga - rUaAf*SfGradp/rhoa);
            Ua.correctBoundaryConditions();

            Ub = HbyAb + fvc::reconstruct(phiDragb - rUbAf*SfGradp/rhob);
            Ub.correctBoundaryConditions();

            U = alpha*Ua + beta*Ub;
        }
    }
}

#include "continuityErrs.H"
