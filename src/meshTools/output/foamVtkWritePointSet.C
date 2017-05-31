/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "foamVtkWritePointSet.H"
#include "foamVtkOutputOptions.H"
#include "OFstream.H"
#include "primitiveMesh.H"
#include "pointSet.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::foamVtkOutput::writePointSet
(
    const primitiveMesh& mesh,
    const pointSet& set,
    const fileName& baseName,
    const foamVtkOutput::outputOptions outOpts
)
{
    outputOptions opts(outOpts);
    opts.legacy(true);  // Legacy only, no append

    std::ofstream os((baseName + (opts.legacy() ? ".vtk" : ".vtp")).c_str());

    autoPtr<foamVtkOutput::formatter> format = opts.newFormatter(os);

    if (opts.legacy())
    {
        foamVtkOutput::legacy::fileHeader(format(), set.name())
            << "DATASET POLYDATA" << nl;
    }

    //-------------------------------------------------------------------------

    const labelList pointLabels(set.sortedToc());

    // Write points
    os << "POINTS " << pointLabels.size() << " float" << nl;

    foamVtkOutput::writeList(format(), mesh.points(), pointLabels);
    format().flush();


    // Write data - pointID
    foamVtkOutput::legacy::pointDataHeader(os, pointLabels.size(), 1);

    os << "pointID 1 " << pointLabels.size() << " int" << nl;

    foamVtkOutput::writeList(format(), pointLabels);
    format().flush();
}


// ************************************************************************* //
