mrfZones.correctBoundaryVelocity(U1);
mrfZones.correctBoundaryVelocity(U2);
mrfZones.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, U2.dimensions()*dimVol/dimTime);

volScalarField dragCoeff(fluid.dragCoeff());

{
    volVectorField liftForce(fluid.liftForce(U));

    {
        U1Eqn =
        (
            fvm::ddt(alpha1, U1)
          + fvm::div(alphaPhi1, U1)

            // Compressibity correction
          - fvm::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), U1)

          + fluid.Cvm()*rho2*alpha1*alpha2/rho1*
            (
                fvm::ddt(U1)
              + fvm::div(phi1, U1)
              - fvm::Sp(fvc::div(phi1), U1)
            )

          + turbulence1->divDevReff(U1)
         ==
          - fvm::Sp(dragCoeff/rho1, U1)
          - alpha1*alpha2/rho1*(liftForce - fluid.Cvm()*rho2*DDtU2)
        );
        mrfZones.addCoriolis(alpha1*(1 + fluid.Cvm()*rho2*alpha2/rho1), U1Eqn);
        U1Eqn.relax();
    }

    {
        U2Eqn =
        (
            fvm::ddt(alpha2, U2)
          + fvm::div(alphaPhi2, U2)

            // Compressibity correction
          - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), U2)

          + fluid.Cvm()*rho2*alpha1*alpha2/rho2*
            (
                fvm::ddt(U2)
              + fvm::div(phi2, U2)
              - fvm::Sp(fvc::div(phi2), U2)
            )
          + turbulence2->divDevReff(U2)
         ==
          - fvm::Sp(dragCoeff/rho2, U2)
          + alpha1*alpha2/rho2*(liftForce + fluid.Cvm()*rho2*DDtU1)
        );
        mrfZones.addCoriolis(alpha2*(1 + fluid.Cvm()*rho2*alpha1/rho2), U2Eqn);
        U2Eqn.relax();
    }
}
