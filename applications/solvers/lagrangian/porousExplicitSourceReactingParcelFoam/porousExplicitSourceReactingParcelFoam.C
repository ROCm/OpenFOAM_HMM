/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    trackedReactingParcelFoam

Description
    - reacting parcel cloud tracking
    - porous media
    - point mass sources
    - polynomial based, incompressible thermodynamics (f(T))

    Note: ddtPhiCorr not used here - not well defined for porous calcs

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hReactionThermo.H"
#include "turbulenceModel.H"
#include "BasicReactingCloud.H"
#include "rhoChemistryModel.H"
#include "chemistrySolver.H"
#include "thermoPhysicsTypes.H"
#include "radiationModel.H"
#include "porousZones.H"
#include "timeActivatedExplicitMulticomponentPointSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readChemistryProperties.H"
    #include "readEnvironmentalProperties.H"
    #include "createFields.H"
    #include "createRadiationModel.H"
    #include "createClouds.H"
    #include "createMulticomponentPointSources.H"
    #include "createPorousZones.H"
    #include "readPISOControls.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readPISOControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        parcels.info();

        #include "chemistry.H"
        #include "rhoEqn.H"
        #include "UEqn.H"

        // --- PIMPLE loop
        for (int oCorr=1; oCorr<=nOuterCorr; oCorr++)
        {
            #include "YEqn.H"
            #include "hEqn.H"

            // --- PISO loop
            for (int corr=1; corr<=nCorr; corr++)
            {
                #include "pEqn.H"
            }

            Info<< "T gas min/max   = " << min(T).value() << ", "
                << max(T).value() << endl;
        }

        turbulence->correct();

        rho = thermo.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
