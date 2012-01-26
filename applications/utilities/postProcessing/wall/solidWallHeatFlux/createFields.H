autoPtr<basicSolidThermo> thermo
(
    basicSolidThermo::New(mesh)
);

volScalarField& T = thermo->T();
