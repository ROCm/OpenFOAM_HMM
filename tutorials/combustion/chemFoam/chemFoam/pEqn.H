{
    thermo.correct();
    rho = thermo.rho();
    if (constProp == "volume")
    {
        p[0] = rho0*R0*thermo.T()[0];
        rho[0] = rho0;
    }
}