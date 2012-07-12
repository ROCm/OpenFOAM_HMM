autoPtr<solidThermo> thermo
(
    solidThermo::New(mesh)
);

const volScalarField& T = thermo->T();
