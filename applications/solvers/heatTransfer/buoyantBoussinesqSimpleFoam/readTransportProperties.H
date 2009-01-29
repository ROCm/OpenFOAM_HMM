    singlePhaseTransportModel laminarTransport(U, phi);

    // thermal expansion coefficient [1/K]
    dimensionedScalar beta(laminarTransport.lookup("beta"));

    // reference temperature [K]
    dimensionedScalar TRef(laminarTransport.lookup("TRef"));

    // reference kinematic pressure [m2/s2]
    dimensionedScalar pRef(laminarTransport.lookup("pRef"));

    // Prandtl number
    dimensionedScalar Pr(laminarTransport.lookup("Pr"));
