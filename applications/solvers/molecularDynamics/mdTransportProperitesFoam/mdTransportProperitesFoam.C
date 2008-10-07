/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2005 OpenCFD Ltd.
    \\/      M anipulation   |
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
    mdTransportProperitesFoam

Description
    MD simulation to calculate continuum transport properites of a homogeneous
    fluid at a given, stationary state.  Density and temperature defined
    by preprocessing, pressure measured.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvCFD.H"
#include "md.H"

int main(int argc, char *argv[])
{
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    moleculeCloud molecules(mesh);

    molecules.removeHighEnergyOverlaps();

#   include "temperatureAndPressureVariables.H"

#   include "createAutoCorrelationFunctions.H"

    label nAveragingSteps = 0;

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime++;

        nAveragingSteps++;

        Info << "Time = " << runTime.timeName() << endl;

        molecules.integrateEquationsOfMotion();

#       include "meanMomentumEnergyAndNMols.H"

#       include "temperatureAndPressure.H"

#       include "calculateAutoCorrelationFunctions.H"

        runTime.write();

        if (runTime.outputTime())
        {
            nAveragingSteps = 0;
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

#   include "calculateTransportProperties.H"

    Info << "End\n" << endl;

    return(0);
}
