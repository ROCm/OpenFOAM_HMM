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
    tmp<volScalarField> tinvTau0
    (
        new volScalarField
        (
            IOobject
            (
                "invTau0",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            invTau
        )
    );
    const volScalarField& invTau0 = tinvTau0();


    // Flow time scale
    // ~~~~~~~~~~~~~~~
    {
        invTau =
            fvc::surfaceSum
            (
                mag(phi)*mesh.deltaCoeffs()/(maxCo*mesh.magSf())
            )
           /rho;

        invTau.max(1.0/maxDeltaT);

        Info<< "    Flow        = "
            << gMin(1/invTau.internalField()) << ", "
            << gMax(1/invTau.internalField()) << endl;
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

        invTau.internalField() = max(invTau.internalField(), 1/tau);
    }


    // Limit rate of change of time scale
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // - reduce as much as required for flow, but limit source contributions
    const dimensionedScalar deltaTRamp("deltaTRamp", dimless, 1/(1 + 0.2));
    invTau = max(invTau, invTau0*deltaTRamp);
    tinvTau0.clear();


    // Limit the largest time scale
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    invTau.max(1/maxDeltaT);


    // Spatially smooth the time scale field
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fvc::smooth(invTau, alphaTauSmooth);

    Info<< "    Overall     = " << min(1/invTau).value()
        << ", " << max(1/invTau).value() << nl << endl;
}


// ************************************************************************* //
