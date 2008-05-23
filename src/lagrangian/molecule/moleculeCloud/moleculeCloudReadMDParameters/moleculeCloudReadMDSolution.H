Info<< nl <<  "Reading MD solution parameters:" << endl;

IOdictionary mdSolution
(
    IOobject
    (
        "mdSolution",
        mesh_.time().system(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

integrationMethod_ = integrationMethodNames_.read
(
    mdSolution.lookup("integrationMethod")
);
Info<< "integrationMethod =\t\t\t"
    << integrationMethodNames_[integrationMethod_] << endl;

potentialEnergyLimit_ = readScalar(mdSolution.lookup("potentialEnergyLimit"));
Info<< "potentialEnergyLimit =\t\t\t" << potentialEnergyLimit_ << endl;

guardRadius_ = readScalar(mdSolution.lookup("guardRadius"));
Info<< "guardRadius =\t\t\t\t" << guardRadius_ << endl;

