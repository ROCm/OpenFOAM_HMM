{
    tmp<volScalarField> tcp(thermo->Cp());
    const volScalarField& cp = tcp();
    rhoCp = rho*cp;

    kappaEff = thermo->kappa() + rho*cp*turbulence->nut()/Prt;

    pDivU = dimensionedScalar("pDivU", p.dimensions()/dimTime, Zero);

DebugVar(fvc::interpolate(rho));

    if (thermo->pDivU())
    {
        pDivU = (p*fvc::div(rhoPhi/fvc::interpolate(rho)));
    }

    const surfaceScalarField rhoCpPhi(fvc::interpolate(cp)*rhoPhi);

    Pair<tmp<volScalarField>> vDotAlphal = mixture->mDot();

    const volScalarField& vDotcAlphal = vDotAlphal[0]();
    const volScalarField& vDotvAlphal = vDotAlphal[1]();
    const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);

    fvScalarMatrix TEqn
    (
         fvm::ddt(rhoCp, T)
       + fvm::div(rhoCpPhi, T, "div(rhoCpPhi,T)")
       - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
       - fvm::laplacian(kappaEff, T)
       + thermo->hc()*vDotvmcAlphal
       + pDivU
    );

    TEqn.relax();
    TEqn.solve();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
