/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    chtMultiRegionFoam

Description
    Combination of heatConductionFoam and buoyantFoam for conjugate heat
    transfer between a solid region and fluid region

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicRhoThermo.H"
#include "turbulenceModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "basicSolidThermo.H"
#include "radiationModel.H"
#include "fvFieldReconstructor.H"
#include "mixedFvPatchFields.H"
#include "fvFieldDecomposer.H"
#include "harmonic.H"
#include "rmap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    regionProperties rp(runTime);

    const label nAllRegions =
        rp.fluidRegionNames().size()
      + rp.solidRegionNames().size();
    PtrList<labelIOList> cellProcAddressing(nAllRegions);
    PtrList<labelIOList> faceProcAddressing(nAllRegions);
    PtrList<labelIOList> boundaryProcAddressing(nAllRegions);
    PtrList<fvMesh> procMeshes(nAllRegions);

    // Load meshes, fluid first
    labelList fluidToProc(identity(rp.fluidRegionNames().size()));
    labelList solidToProc(rp.solidRegionNames().size());
    forAll(solidToProc, i)
    {
        solidToProc[i] = fluidToProc.size()+i;
    }

    // Get the coupled solution flag
    #include "readPIMPLEControls.H"

    if (temperatureCoupled)
    {
        Info<< "Solving single enthalpy for all equations" << nl << endl;
    }

    #include "createFluidMeshes.H"
    #include "createSolidMeshes.H"

    #include "createFluidFields.H"
    #include "createSolidFields.H"

    // Temperature solved on single mesh
    #include "createAllMesh.H"
    #include "createAllFields.H"

    #include "initContinuityErrs.H"

    #include "readTimeControls.H"
    #include "readSolidTimeControls.H"


    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"


    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"


        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;


        if (nOuterCorr != 1)
        {
            forAll(fluidToProc, i)
            {
                #include "setRegionFluidFields.H"
                #include "storeOldFluidFields.H"
            }
        }


        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        {
            bool finalIter = oCorr == nOuterCorr-1;

            if (finalIter)
            {
                forAll(procMeshes, procI)
                {
                    procMeshes[procI].data::add("finalIteration", true);
                }
            }


            PtrList<surfaceScalarField> procPhi(nAllRegions);
            PtrList<surfaceScalarField> procAlpha(nAllRegions);


            // Solve (uncoupled) or set up (coupled) the temperature equation
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            forAll(solidToProc, i)
            {
                label procI = solidToProc[i];

                Info<< "\nSolving temperature for solid region "
                    << procMeshes[procI].name() << endl;
                #include "setRegionSolidFields.H"
                #include "readSolidMultiRegionPIMPLEControls.H"

                if (temperatureCoupled)
                {
                    // Map my properties to overall h equation
                    #include "rmapSolid.H"
                }
                else
                {
                    #include "solveSolid.H"
                }
            }


            forAll(fluidToProc, i)
            {
                label procI = fluidToProc[i];

                Info<< "\nSolving temperature for fluid region "
                    << procMeshes[procI].name() << endl;
                #include "setRegionFluidFields.H"
                #include "readFluidMultiRegionPIMPLEControls.H"

                if (oCorr == 0)
                {
                    #include "rhoEqn.H"
                }

                if (temperatureCoupled)
                {
                    // Map my properties to overall h equation
                    #include "rmapFluid.H"
                }
                else
                {
                    #include "hEqn.H"
                }
            }


            // Solve combined h equation
            // ~~~~~~~~~~~~~~~~~~~~~~~~~

            if (temperatureCoupled)
            {
                Info<< "\nSolving single enthalpy for all regions"
                    << endl;

                // Solve combined h
                #include "allhEqn.H"

                forAll(solidToProc, i)
                {
                    label procI = solidToProc[i];
                    #include "setRegionSolidFields.H"
                    #include "readSolidMultiRegionPIMPLEControls.H"

                    #include "mapSolid.H"
                }

                forAll(fluidToProc, i)
                {
                    label procI = fluidToProc[i];
                    #include "setRegionFluidFields.H"
                    #include "readFluidMultiRegionPIMPLEControls.H"

                    #include "mapFluid.H"
                }
            }


            // Update thermos
            // ~~~~~~~~~~~~~~

            forAll(fluidToProc, i)
            {
                Info<< "\nUpdating thermo for fluid region "
                    << procMeshes[fluidToProc[i]].name() << endl;

                #include "setRegionFluidFields.H"
                #include "readFluidMultiRegionPIMPLEControls.H"

                thermo.correct();
                rad.correct();
                #include "solvePressureVelocityFluid.H"
            }

            forAll(solidToProc, i)
            {
                Info<< "\nUpdating thermo for solid region "
                    << procMeshes[solidToProc[i]].name() << endl;
                #include "setRegionSolidFields.H"
                #include "readSolidMultiRegionPIMPLEControls.H"

                thermo.correct();
            }


            if (finalIter)
            {
                forAll(procMeshes, procI)
                {
                    procMeshes[procI].data::remove("finalIteration");
                }
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
