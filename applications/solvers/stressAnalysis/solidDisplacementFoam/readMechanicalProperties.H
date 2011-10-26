    Info<< "Reading mechanical properties\n" << endl;

    IOdictionary mechanicalProperties
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const dictionary& rhoDict(mechanicalProperties.subDict("rho"));
    word rhoType(rhoDict.lookup("rho"));

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimVolume, 0.0)
    );

    if (rhoType == "rhoInf")
    {
        rho = rhoDict.lookup("rhoInf");
    }

    volScalarField rhoE
    (
        IOobject
        (
            "E",
            runTime.timeName(0),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimMass/dimLength/sqr(dimTime), 0.0)
    );

    const dictionary& EDict(mechanicalProperties.subDict("E"));
    word EType(EDict.lookup("E"));
    if (EType == "EInf")
    {
        rhoE = EDict.lookup("EInf");
    }


    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(0),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    );

    const dictionary& nuDict(mechanicalProperties.subDict("nu"));
    word nuType(nuDict.lookup("nu"));

    if (nuType == "nuInf")
    {
        nu = nuDict.lookup("nuInf");
    }

    Info<< "Normalising E : E/rho\n" << endl;
    volScalarField E = rhoE/rho;

    Info<< "Calculating Lame's coefficients\n" << endl;

    volScalarField mu = E/(2.0*(1.0 + nu));
    volScalarField lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));
    volScalarField threeK = E/(1.0 - 2.0*nu);

    Switch planeStress(mechanicalProperties.lookup("planeStress"));

    if (planeStress)
    {
        Info<< "Plane Stress\n" << endl;

        lambda = nu*E/((1.0 + nu)*(1.0 - nu));
        threeK = E/(1.0 - nu);
    }
    else
    {
        Info<< "Plane Strain\n" << endl;
    }
