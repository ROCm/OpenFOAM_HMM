{
    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phi1 - phi2);

    surfaceScalarField alpha1f(fvc::interpolate(max(alpha1, scalar(0))));

    tmp<surfaceScalarField> pPrimeByA;

    if (implicitPhasePressure)
    {
        pPrimeByA =
            fvc::interpolate((1.0/rho1)*rAU1*turbulence1().pPrime())
          + fvc::interpolate((1.0/rho2)*rAU2*turbulence2().pPrime());

        surfaceScalarField phiP
        (
            pPrimeByA()*fvc::snGrad(alpha1, "bounded")*mesh.magSf()
        );

        phic += alpha1f*phiP;
        phir += phiP;
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::DimensionedInternalField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("Sp", dgdt.dimensions(), 0.0)
        );

        volScalarField::DimensionedInternalField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi)*min(alpha1, scalar(1))
        );

        forAll(dgdt, celli)
        {
            if (dgdt[celli] > 0.0 && alpha1[celli] > 0.0)
            {
                Sp[celli] -= dgdt[celli]*alpha1[celli];
                Su[celli] += dgdt[celli]*alpha1[celli];
            }
            else if (dgdt[celli] < 0.0 && alpha1[celli] < 1.0)
            {
                Sp[celli] += dgdt[celli]*(1.0 - alpha1[celli]);
            }
        }

        dimensionedScalar totalDeltaT = runTime.deltaT();
        if (nAlphaSubCycles > 1)
        {
            alphaPhi1 = dimensionedScalar("0", alphaPhi1.dimensions(), 0);
        }

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            surfaceScalarField alphaPhic1
            (
                fvc::flux
                (
                    phic,
                    alpha1,
                    alphaScheme
                )
              + fvc::flux
                (
                    -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                    alpha1,
                    alpharScheme
                )
            );

            alphaPhic1.boundaryField() = min
            (
                phi1.boundaryField()*alpha1.boundaryField(),
                alphaPhic1.boundaryField()
            );

            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                alphaPhic1,
                Sp,
                Su,
                1,
                0
            );

            if (nAlphaSubCycles > 1)
            {
                alphaPhi1 += (runTime.deltaT()/totalDeltaT)*alphaPhic1;
            }
            else
            {
                alphaPhi1 = alphaPhic1;
            }
        }

        if (implicitPhasePressure)
        {
            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1) - fvc::ddt(alpha1)
              - fvm::laplacian(alpha1f*pPrimeByA, alpha1, "bounded")
            );

            alpha1Eqn.relax();
            alpha1Eqn.solve();

            alphaPhi1 += alpha1Eqn.flux();
        }

        alphaPhi2 = phi - alphaPhi1;
        alpha2 = scalar(1) - alpha1;

        Info<< "Dispersed phase volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;
    }
}

rho = fluid.rho();
