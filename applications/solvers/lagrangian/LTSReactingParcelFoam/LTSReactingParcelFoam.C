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

Application
    LTSReactingParcelFoam

Description
    Local time stepping (LTS) solver for steady, compressible, laminar or
    turbulent reacting and non-reacting flow with multiphase Lagrangian
    parcels and porous media, including explicit sources for mass, momentum
    and energy

    Note: ddtPhiCorr not used here when porous zones are active
    - not well defined for porous calculations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hReactionThermo.H"
#include "turbulenceModel.H"
#include "basicReactingMultiphaseCloud.H"
#include "rhoChemistryModel.H"
#include "chemistrySolver.H"
#include "radiationModel.H"
#include "porousZones.H"
#include "timeActivatedExplicitSource.H"
#include "SLGThermo.H"
#include "fvcSmooth.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"

    pimpleControl pimple(mesh);

    #include "readTimeControls.H"
    #include "readAdditionalSolutionControls.H"
    #include "createFields.H"
    #include "createRadiationModel.H"
    #include "createClouds.H"
    #include "createExplicitSources.H"
    #include "createPorousZones.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readChemistryProperties.H"
        #include "readAdditionalSolutionControls.H"
        #include "readTimeControls.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        #include "chemistry.H"
        #include "timeScales.H"

        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        for (pimple.start(); pimple.loop(); pimple++)
        {
            if (pimple.nOuterCorr() != 1)
            {
                p.storePrevIter();
            }

            turbulence->correct();

            #include "UEqn.H"
            #include "YEqn.H"
            #include "hsEqn.H"

            // --- PISO loop
            for (int corr=0; corr<pimple.nCorr(); corr++)
            {
                #include "pEqn.H"
            }
        }

        if (runTime.write())
        {
            chemistry.dQ()().write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
