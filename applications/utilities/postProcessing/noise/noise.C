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

Description
    Utility to perform noise analysis of pressure data.


    \heading Usage

    Control settings are read from the $FOAM_CASE/system/noiseDict dictionary,
    or user-specified dictionary using the -dict option.  Presure data can be
    supplied as a single time history at a given point location, or a surface
    of time histories, based on the run-time selectable \c noiseModel, e.g.

    \vebatim
    noiseModel      <geometryType>Noise;

    <geometryType>Coeffs
    {
        ...
    }
    \endverbatim

    where \<geometryMode\> is either \c point or \surface

    Both operation types require:
    \verbatim
    pRef            101325;
    N               65536;
    fl              25;
    fu              10000;
    fftWriteInterval 1;
    dataMode        point;
    \endverbatim

    where
    \table
        Property    | Description                   | Required  | Default value
        pRef        | Reference pressure            | no        | 0
        N           | Number of samples in sampling window | no | 65536 (2^16)
        fl          | Lower frequency band          | no        | 25
        fu          | Upper frequency band          | no        | 10000
        fftWriteInterval | Output interval for FFT data | no    | 1
        dataMode    | Input data mode, 'point' or 'surface' | yes |
    \endtable

    Point-based pressure data is read using a CSV reader, and output in a
    graph format:
    \verbatim
    pointData
    {
        csvFileData
        {
            fileName        "pressureData";
            nHeaderLine     1;
            refColumn       0;
            componentColumns (1);
            separator       " ";
            mergeSeparators yes;
        }

        graphFormat     raw;
        windowOverlapPercent 0;
        nWindow         100; // fixed number of windows
    }
    \endverbatim

    where
    \table
        Property    | Description                   | Required  | Default value
        csvFileData | CSV file data dictionary - see CSV.H  | yes |
        graphFormat | Output graph format           | no        | raw
        windowOverlapPercent | Window overla percent| yes       |
        nWindow     | Number of sampling windows    | no        | 1, calculated
    \endtable

    Surface-based pressure data is specified using a surface reader, and output
    using a surface writer.  The reader type must support multiple time
    handling, e.g. the ensight format.  Output surfaces are written on a
    per-frequency basis:
    \verbatim
    surfaceData
    {
        reader          ensight;
        writer          ensight;
        inputFile       "postProcessing/faceSource1/surface/patch1/patch1.case";
        pName           p;
        windowOverlapPercent 50;
        // nWindow         0; // let code decide if not present
    }
    \endverbatim

    where
    \table
        Property    | Description                   | Required  | Default value
        reader      | Surface reader type           | yes       |
        writer      | Surface writer type           | yes       |
        inputFile   | Pressure data surface file    | yes       |
        pName       | Name of pressure field        | no        | p
        windowOverlapPercent | Window overla percent| yes       |
        nWindow     | Number of sampling windows    | no        | 1, calculated
    \endtable

    Current outputs include:
    - FFT of the pressure data (Pf)
    - narrow-band PFL (pressure-fluctuation level) spectrum (Lf)
    - one-third-octave-band PFL spectrum (Ldelta)
    - one-third-octave-band pressure spectrum (Pdelta)

Note:
  - If using the windowOverlapPercent entry, the code will automatically
    determine number of windows
  - However, if supplied, the nWindow entry will take precedence

SeeAlso
    CSV.H
    surfaceReader.H
    surfaceWriter.H
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
