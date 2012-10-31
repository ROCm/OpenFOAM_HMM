    MRFZones mrfZones(mesh);
    mrfZones.correctBoundaryVelocity(U);

    IOporosityModelList pZones(mesh);
    Switch pressureImplicitPorosity(false);
