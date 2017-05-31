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

#include "foamVtkWriteCellSetFaces.H"
#include "foamVtkOutputOptions.H"
#include "OFstream.H"
#include "primitiveMesh.H"
#include "cellSet.H"
#include "uindirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::foamVtkOutput::writeCellSetFaces
(
    const primitiveMesh& mesh,
    const cellSet& set,
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

    // External faces of cellset with OpenFOAM cellID as value

    Map<label> cellFaces(2*set.size());

    forAllConstIters(set, iter)
    {
        label celli = iter.key();
        const cell& cFaces = mesh.cells()[celli];

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if (mesh.isInternalFace(facei))
            {
                label otherCelli = mesh.faceOwner()[facei];

                if (otherCelli == celli)
                {
                    otherCelli = mesh.faceNeighbour()[facei];
                }

                if (!set.found(otherCelli))
                {
                    cellFaces.insert(facei, celli);
                }
            }
            else
            {
                cellFaces.insert(facei, celli);
            }
        }
    }

    const labelList faceLabels = cellFaces.sortedToc();
    labelList faceValues(cellFaces.size());

    forAll(faceLabels, facei)
    {
        faceValues[facei] = cellFaces[faceLabels[facei]];  // Cell ID
    }

    uindirectPrimitivePatch pp
    (
        UIndirectList<face>(mesh.faces(), faceLabels),
        mesh.points()
    );

    //-------------------------------------------------------------------------

    // Write points and faces as polygons
    os << "POINTS " << pp.nPoints() << " float" << nl;

    foamVtkOutput::writeList(format(), pp.localPoints());
    format().flush();

    label count = pp.size();
    forAll(pp, facei)
    {
        count += pp.localFaces()[facei].size();
    }
    os << "POLYGONS " << pp.size() << ' ' << count << nl;

    forAll(pp, facei)
    {
        const face& f = pp.localFaces()[facei];

        format().write(f.size());
        foamVtkOutput::writeList(format(), f);
    }
    format().flush();


    // Write data - faceId/cellId
    foamVtkOutput::legacy::cellDataHeader(os, pp.size(), 1);

    os << "cellID 1 " << pp.size() << " int" << nl;

    foamVtkOutput::writeList(format(), faceValues);
    format().flush();
}


// ************************************************************************* //
