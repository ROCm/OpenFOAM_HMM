Info<< "\nConstructing interpolation" << endl;

Info << "\nConstructing gas properties" << endl;
/*
PtrList<specieConstProperties> gasProperties(Y.size());
forAll(gasProperties, i)
{
    gasProperties.set
    (
        i,
        new specieConstProperties
        (
            dynamic_cast<const multiComponentMixture<constTransport<
                specieThermo<hConstThermo<perfectGas> > > >&>
                (thermo()).speciesData()[i]
        )
    );
}
*/
PtrList<gasThermoPhysics> gasProperties(Y.size());
forAll(gasProperties, i)
{
    gasProperties.set
    (
        i,
        new gasThermoPhysics
        (
            dynamic_cast<const reactingMixture<gasThermoPhysics>&>(thermo).speciesData()[i]
        )
    );
}

Info<< "\nConstructing coal cloud" << endl;
coalCloud coalParcels
(
    "coalCloud1",
    rho,
    U,
    g,
    thermo,
    gasProperties
);

Info<< "\nConstructing limestone cloud" << endl;
basicThermoCloud limestoneParcels
(
    "limestoneCloud1",
    rho,
    U,
    g,
    thermo
);
