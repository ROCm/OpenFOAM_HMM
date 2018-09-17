Info<< "\nConstructing " << CLOUD_BASE_TYPE_NAME << " cloud" << endl;
basicReactingTypeCloud parcels
(
    word(CLOUD_BASE_TYPE_NAME) + "1",
    rho,
    U,
    g,
    slgThermo
);
