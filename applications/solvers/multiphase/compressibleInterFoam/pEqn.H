{
    rho1 = rho10 + psi1*p;
    rho2 = rho20 + psi2*p;

    volScalarField rAU = 1.0/UEqn.A();
    surfaceScalarField rAUf = fvc::interpolate(rAU);

    tmp<fvScalarMatrix> p_rghEqnComp1;
    tmp<fvScalarMatrix> p_rghEqnComp2;

    //if (transonic)
    //{
    //}
    //else
    {
        surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi);
        surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi);

        p_rghEqnComp1 =
            fvc::ddt(rho1) + psi1*correction(fvm::ddt(p_rgh))
          + fvc::div(phid1, p_rgh)
          - fvc::Sp(fvc::div(phid1), p_rgh);

        p_rghEqnComp2 =
            fvc::ddt(rho2) + psi2*correction(fvm::ddt(p_rgh))
          + fvc::div(phid2, p_rgh)
          - fvc::Sp(fvc::div(phid2), p_rgh);
    }

    U = rAU*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

    phi = phiU +
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    //thermo.rho() -= psi*p_rgh;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phi)
          - fvm::laplacian(rAUf, p_rgh)
        );

        solve
        (
            (
                (max(alpha1, scalar(0))/rho1)*p_rghEqnComp1()
              + (max(alpha2, scalar(0))/rho2)*p_rghEqnComp2()
            )
          + p_rghEqnIncomp,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            // Second part of thermodynamic density update
            //thermo.rho() += psi*p_rgh;

            dgdt =
            (
                pos(alpha2)*(p_rghEqnComp2 & p_rgh)/rho2
              - pos(alpha1)*(p_rghEqnComp1 & p_rgh)/rho1
            );
            phi += p_rghEqnIncomp.flux();
        }
    }

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();

    p = max
    (
        (p_rgh + gh*(alpha1*rho10 + alpha2*rho20))
       /(1.0 - gh*(alpha1*psi1 + alpha2*psi2)),
        pMin
    );

    rho1 = rho10 + psi1*p;
    rho2 = rho20 + psi2*p;

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;
}
