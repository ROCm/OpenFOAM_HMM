/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "volPointInterpolation.H"
#include "PtrList.H"

#include "writer.H"
#include "sampleSet.H"
#include "volFieldSampler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"


    //
    // 1. Sampling and writing.
    //

    // Set up mesh searching
    meshSearch searchEngine(mesh, true);

    // Construct sample points generator
    dictionary sampleDict;
    sampleDict.add("name", "lineX1");
    sampleDict.add("axis", "distance");
    sampleDict.add("start", point(0.02, 0.051, 0.005));
    sampleDict.add("end", point(0.06, 0.051, 0.005));
    sampleDict.add("nPoints", 10);

    PtrList<sampleSet> sampleSets(1);
    sampleSets.set
    (
        0,
        sampleSet::New
        (
            "uniform",
            mesh,
            searchEngine,
            sampleDict
        ).ptr()
    );

    // Load field
    volScalarField sField
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );

    // Set up interpolation
    pointMesh pMesh(mesh);
    volPointInterpolation pInterp(mesh, pMesh);


    // Do actual sampling
    dictionary interpolationSchemes;
    interpolationSchemes.add("p", "cell");

    volFieldSampler<scalar> sampleVals
    (
        pInterp,
        interpolationSchemes,
        sField,
        sampleSets
    );

    // Construct writer and write
    autoPtr<writer<scalar> > scalarFormatter(writer<scalar>::New("xmgr"));

    scalarFormatter().write
    (
        ".",
        sampleSets[0],
        sampleVals.name(),
        sampleVals[0]
    );



    //
    // 2. No sampling, just using writing
    //
    List<point> points(5);
    points[0] = point(0, 0, 0);
    points[1] = point(1, 1, 1);
    points[2] = point(2, 2, 2);
    points[3] = point(3, 3, 3);
    points[4] = point(4, 4, 4);

    scalarList vals(5);
    vals[0] = 0.0;
    vals[1] = 0.1;
    vals[2] = 0.2;
    vals[3] = 0.3;
    vals[4] = 0.4;


    scalarFormatter().write
    (
        ".",
        coordSet
        (
            "someLine",
            "distance",
            points,
            points[0]
        ),
        "U.component(0)",
        vals
    );


    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
