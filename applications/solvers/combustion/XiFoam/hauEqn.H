if (ign.ignited())
{
    solve
    (
        fvm::ddt(rho, hau)
      + mvConvection->fvmDiv(phi, hau)
      - fvm::laplacian(turbulence->alphaEff(), hau)

    // These terms cannot be used in partially-premixed combustion due to
    // the resultant inconsistency between ft and hau transport.
    // A possible solution would be to solve for ftu as well as ft.
    //- fvm::div(muEff*fvc::grad(b)/(b + 0.001), hau)
    //+ fvm::Sp(fvc::div(muEff*fvc::grad(b)/(b + 0.001)), hau)

    ==
        (dpdt - (fvc::ddt(rho, K) + fvc::div(phi, K)))*rho/thermo.rhou()
    );
}
