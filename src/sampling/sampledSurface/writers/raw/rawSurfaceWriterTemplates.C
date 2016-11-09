/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "OFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::rawSurfaceWriter::writeTemplate
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
    const pointField& points = surf.points();
    const faceList&    faces = surf.faces();

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os(outputDir/fieldName + '_' + surfaceName + ".raw");

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << os.name() << endl;
    }

    // header
    os  << "# " << fieldName;
    if (isNodeValues)
    {
        os  << "  POINT_DATA ";
    }
    else
    {
        os  << "  FACE_DATA ";
    }

    // header
    writeHeader(os, fieldName, values);

    // values
    if (isNodeValues)
    {
        forAll(values, elemI)
        {
            writeLocation(os, points[elemI]);
            writeData(os, values[elemI]);
        }
    }
    else
    {
        forAll(values, elemI)
        {
            writeLocation(os, points, faces[elemI]);
            writeData(os, values[elemI]);
        }
    }

    return os.name();
}


// ************************************************************************* //
