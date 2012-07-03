autoPtr<solidThermo> thermo
(
    solidThermo::New(mesh)
);

volScalarField& T = thermo->T();
