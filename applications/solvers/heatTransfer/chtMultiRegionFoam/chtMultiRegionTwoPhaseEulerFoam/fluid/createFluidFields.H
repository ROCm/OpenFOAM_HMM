// Initialise fluid field pointer lists
PtrList<twoPhaseSystem> phaseSystemFluid(fluidRegions.size());

PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());

PtrList<volScalarField> p_rghFluid(fluidRegions.size());

PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>
    fieldsFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());
List<bool> frozenFlowFluid(fluidRegions.size(), false);

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());
PtrList<dimensionedScalar> pMinFluid(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

PtrList<pimpleControl> pimpleFluid(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    pimpleFluid.set
    (
        i,
        new pimpleControl(fluidRegions[i])
    );

    Info<< "    Adding to phaseSystemFluid\n" << endl;
    phaseSystemFluid.set(i, twoPhaseSystem::New(fluidRegions[i]).ptr());

    Info<< "    Adding hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, Zero)
        )
    );

    Info<< "    Adding ghRef\n" << endl;
    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding p_rghFluid\n" << endl;
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Correcting p_rghFluid\n" << endl;

    // Force p_rgh to be consistent with p
    p_rghFluid[i] =
        phaseSystemFluid[i].phase1().thermo().p()
     -  phaseSystemFluid[i].phase1().thermo().rho()*ghFluid[i];

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    Info<< "    Correcting initialMassFluid\n" << endl;
    initialMassFluid[i] =
        fvc::domainIntegrate(phaseSystemFluid[i].rho()).value();

    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");

    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    pRefCellFluid[i] = -1;
    pRefValueFluid[i] = 0.0;

    Info<< "    Setting reference\n" << endl;
    if (p_rghFluid[i].needReference())
    {
        setRefCell
        (
            phaseSystemFluid[i].phase1().thermoRef().p(),
            p_rghFluid[i],
            pimpleDict,
            pRefCellFluid[i],
            pRefValueFluid[i]
        );
    }

    pMinFluid.set
    (
        i,
        new dimensionedScalar
        (
            "pMin",
            dimPressure,
            phaseSystemFluid[i]
        )
    );
}
