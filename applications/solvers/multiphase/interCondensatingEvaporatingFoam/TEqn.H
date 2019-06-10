{
    tmp<volScalarField> tcp(thermo->Cp());
    const volScalarField& cp = tcp();

    rhoCp = rho*cp;

    kappaEff = thermo->kappa() + rho*cp*turbulence->nut()/Prt;

    const surfaceScalarField rhoCpPhi(fvc::interpolate(cp)*rhoPhi);

    Pair<tmp<volScalarField>> vDotAlphal = mixture->mDot();

    const volScalarField& vDotcAlphal = vDotAlphal[0]();
    const volScalarField& vDotvAlphal = vDotAlphal[1]();
    const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);

    fvScalarMatrix TEqn
    (
         fvm::ddt(rhoCp, T)
       + fvm::div(rhoCpPhi, T)
       - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
       - fvm::laplacian(kappaEff, T)
       + thermo->hc()*vDotvmcAlphal
    );


    TEqn.relax();
    TEqn.solve();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
