    Info<< "Creating radiation model\n" << endl;
    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(thermo.T())
    );
