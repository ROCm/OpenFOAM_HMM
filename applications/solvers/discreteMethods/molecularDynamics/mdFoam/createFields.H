    Info<< "\nReading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    word cloudName("polyatomicCloud");

    potential pot
    (
        mesh,
        IOdictionary
        (
            IOobject
            (
                cloudName +  "Properties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    );

    polyatomicCloud molecules(cloudName, mesh, pot);
