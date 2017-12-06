/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "OFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"
#include "Time.H"
#include "pointIOField.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::boundaryDataSurfaceWriter::writeTemplate
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    const fileName baseDir(outputDir.path()/surfaceName);
    const fileName timeName(outputDir.name());

    const pointField& points = surf.points();
    const faceList&    faces = surf.faces();


    // Construct dummy time to use as an objectRegistry
    const fileName caseDir(getEnv("FOAM_CASE"));
    Time dummyTime
    (
        caseDir.path(), //rootPath,
        caseDir.name(), //caseName,
        "system",       //systemName,
        "constant",     //constantName,
        false           //enableFunctionObjects
    );


    // Write points

    pointIOField pts
    (
        IOobject
        (
            baseDir/"points",
            dummyTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        label(0)
    );

    if (isNodeValues)
    {
        if (verbose)
        {
            Info<< "Writing points to " << baseDir/"points" << endl;
        }
        pts = points;
    }
    else
    {
        if (verbose)
        {
            Info<< "Writing face centres to " << baseDir/"points" << endl;
        }

        primitivePatch pp(SubList<face>(faces, faces.size()), points);

        pts = pp.faceCentres();
    }

    {
        // Do like regIOobject::writeObject but don't do instance() adaptation
        // since this would write to e.g. 0/ instead of postProcessing/

        // Try opening an OFstream for object
        mkDir(pts.path());
        OFstream os(pts.objectPath());

        pts.writeHeader(os);
        pts.writeData(os);
        pts.writeEndDivider(os);
    }


    // Write field
    {
        fileName valsDir(baseDir/timeName);
        mkDir(valsDir);
        OFstream os(valsDir/fieldName);
        os  << values;
    }

    return baseDir;
}


// ************************************************************************* //
