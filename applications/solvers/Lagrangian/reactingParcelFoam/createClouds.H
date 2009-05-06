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
PtrList<specieReactingProperties> gasProperties(Y.size());
forAll(gasProperties, i)
{
    gasProperties.set
    (
        i,
        new specieReactingProperties
        (
            dynamic_cast<const reactingMixture&>(thermo()).speciesData()[i]
        )
    );
}

Info<< "\nConstructing reacting cloud" << endl;
basicReactingCloud reactingParcels
(
    "reactingCloud1",
    rho,
    U,
    g,
    thermo(),
    gasProperties
);
