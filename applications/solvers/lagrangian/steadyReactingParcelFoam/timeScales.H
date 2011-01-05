/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

Info<< "Time scales min/max:" << endl;


{
    // Cache old time scale field
    tmp<volScalarField> tinvTauFlow0
    (
        new volScalarField
        (
            IOobject
            (
                "invTauFlow0",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            invTauFlow
        )
    );
    const volScalarField& invTauFlow0 = tinvTauFlow0();


    // Flow time scale
    // ~~~~~~~~~~~~~~~
    {
        invTauFlow =
            fvc::surfaceSum
            (
                mag(phi)*mesh.deltaCoeffs()/(maxCo*mesh.magSf())
            )
           /rho;

        invTauFlow.max(1.0/maxDeltaT);

        Info<< "    Flow        = "
            << gMin(1/invTauFlow.internalField()) << ", "
            << gMax(1/invTauFlow.internalField()) << endl;
    }


    // Mass source time scale
    // ~~~~~~~~~~~~~~~~~~~~~~

    {
        scalarField tau
        (
            runTime.deltaTValue()*mag(parcels.Srho() + massSource.SuTot())
        );

        tau = alphaTauRho*rho/(tau + ROOTVSMALL);

        Info<< "    Density     = "
            << min(maxDeltaT, gMin(tau)) << ", "
            << min(maxDeltaT, gMax(tau)) << endl;

        invTauFlow.internalField() = max(invTauFlow.internalField(), 1/tau);
    }


    // Momentum source time scale
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
/*
        // Method 1 - mag(U) limit using 'small' nominal velocity
        scalarField tau
        (
            runTime.deltaTValue()
           *mag
            (
                rho.dimensionedInternalField()*g
              + parcels.UTrans()/(mesh.V()*runTime.deltaT())
              + momentumSource.Su()
            )
           /rho
        );

        const scalar nomMagU(dimensionedScalar("1", dimVelocity, 1));
        tau = alphaTauU*(nomMagU + mag(U))/(tau + ROOTVSMALL);
*/
/*
        // Method 2 - based on fluxes and Co-like limit
        volVectorField UEqnRhs
        (
            IOobject
            (
                "UEqnRhs",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedVector("zero", dimDensity*dimAcceleration, vector::zero),
            zeroGradientFvPatchVectorField::typeName
        );

        UEqnRhs.internalField() =
            rho.dimensionedInternalField()*g
          + parcels.UTrans()/(mesh.V()*runTime.deltaT())
          + momentumSource.Su();

        surfaceScalarField phiSU
        (
            "phiSU",
            fvc::interpolate(runTime.deltaT()*UEqnRhs) & mesh.Sf()
        );

        scalarField tau
        (
            alphaTauU*rho
           /fvc::surfaceSum
            (
                mag(phi + phiSU)*mesh.deltaCoeffs()/mesh.magSf()
              + dimensionedScalar("SMALL", dimDensity/dimTime, ROOTVSMALL)
            )
        );

*/
/*
        Info<< "    Momentum    = " << min(maxDeltaT, gMin(tau)) << ", "
            << min(maxDeltaT, gMax(tau)) << endl;

        invTauFlow.internalField() = max(invTauFlow.internalField(), 1/tau);
*/
    }


    // Temperature source time scale
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        scalarField tau
        (
            runTime.deltaTValue()
           *mag
            (
                DpDt
              + parcels.hsTrans()/(mesh.V()*runTime.deltaT())
              + energySource.Su()
              + chemistrySh
            )
           /rho
        );

        tau = alphaTauTemp*thermo.Cp()*T/(tau + ROOTVSMALL);

        Info<< "    Temperature = " << min(maxDeltaT, gMin(tau)) << ", "
            << min(maxDeltaT, gMax(tau)) << endl;

        invTauFlow.internalField() = max(invTauFlow.internalField(), 1/tau);
    }


    // Specie source time scale
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    {
        scalarField tau(mesh.nCells(), ROOTVGREAT);
        forAll(Y, fieldI)
        {
            const volScalarField& Yi = Y[fieldI];
            const scalarField deltaYi
            (
                runTime.deltaTValue()
               *mag
                (
                    kappa*chemistry.RR(fieldI)()
                  + massSource.Su(fieldI)
                  + parcels.Srho(fieldI)
                )
               /rho
            );
            tau =
                min
                (
                    tau,
                    alphaTauSpecie
                   /(deltaYi/(Yi + specieMaxUnbound) + ROOTVSMALL)
                );
        }

        Info<< "    Specie      = " << min(maxDeltaT, gMin(tau)) << ", "
            << min(maxDeltaT, gMax(tau)) << endl;

        invTauFlow.internalField() = max(invTauFlow.internalField(), 1/tau);
    }


    // Limit rate of change of time scale
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // - reduce as much as required for flow, but limit source contributions
    const dimensionedScalar deltaTRamp("deltaTRamp", dimless, 1/(1 + 0.2));
    invTauFlow = max(invTauFlow, invTauFlow0*deltaTRamp);
    tinvTauFlow0.clear();


    // Limit the largest time scale
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    invTauFlow.max(1/maxDeltaT);


    // Spatially smooth the time scale field
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fvc::smooth(invTauFlow, alphaTauSmooth);

    Info<< "    Overall     = " << min(1/invTauFlow).value()
        << ", " << max(1/invTauFlow).value() << nl << endl;
}


// ************************************************************************* //
