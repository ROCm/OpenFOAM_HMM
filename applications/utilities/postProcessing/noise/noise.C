/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    noise

Group
    grpPostProcessingUtilities

Description
    Utility to perform noise analysis of pressure data.

    The utility provides a light wrapper around the run-time selectable
    noise model.  Current options include:
    - point, and
    - surface noise.

    \heading Example usage

    \verbatim
    noiseModel      surfaceNoise; // pointNoise

    surfaceNoiseCoeffs
    {
        windowModel     Hanning;

        HanningCoeffs
        {
            // Window overlap percentage
            overlapPercent  50;
            symmetric       yes;
            extended        yes;

            // Optional number of windows, default = all available
            nWindow         5;
        }


        // Input file
        inputFile   "postProcessing/faceSource1/surface/patch/patch.case";

        // Surface reader
        reader      ensight;

        // Surface writer
        writer      ensight;

        // Collate times for ensight output - ensures geometry is only written once
        writeOptions
        {
            ensight
            {
                collateTimes 1;
            }
        }

        // Number of samples in sampling window
        // Must be a power of 2, default = 2^16 (=65536)
        N               4096;

        // Write interval for FFT data, default = 1
        fftWriteInterval 100;
    }
    \endverbatim


SeeAlso
    noiseFFT.H
    noiseModel.H
    windowModel.H

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "noiseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    word dictName("noiseDict");
    if (args.optionFound("dict"))
    {
        dictName = args["dict"];
    }

    IOdictionary dict
    (
        IOobject
        (
            dictName,
            runTime.system(),
            runTime,
            IOobject::MUST_READ
        )
    );

    autoPtr<noiseModel> model(noiseModel::New(dict));
    model->calculate();

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
