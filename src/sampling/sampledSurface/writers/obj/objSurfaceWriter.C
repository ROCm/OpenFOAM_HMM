/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "objSurfaceWriter.H"
#include "fileName.H"
#include "OFstream.H"
#include "faceList.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::objSurfaceWriter<Type>::objSurfaceWriter()
:
    surfaceWriter<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::objSurfaceWriter<Type>::~objSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::objSurfaceWriter<Type>::write
(
    const fileName& samplePath,
    const fileName& timeDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const fileName& fieldName,
    const Field<Type>& values,
    const bool verbose
) const
{
    fileName surfaceDir(samplePath/timeDir);

    if (!isDir(surfaceDir))
    {
        mkDir(surfaceDir);
    }

    fileName fName(surfaceDir/surfaceName + ".obj");

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << fName << endl;
    }

    // this is a quick hack
    OFstream os(fName);

    os  << "# Wavefront OBJ file" << nl
        << "o " << os.name().lessExt().name() << nl
        << nl
        << "# points : " << points.size() << nl
        << "# faces  : " << faces.size() << nl
        << "# no zones " << nl;

    os  << nl
        << "# <points count=\"" << points.size() << "\">" << endl;

    // Write vertex coords
    forAll(points, ptI)
    {
        os  << "v " << points[ptI].x()
            << ' '  << points[ptI].y()
            << ' '  << points[ptI].z() << nl;
    }

    os  << "# </points>" << nl
        << nl
        << "# <faces count=\"" << faces.size() << "\">" << endl;

    forAll(faces, i)
    {
        const face& f = faces[i];

        os << 'f';
        forAll(f, fp)
        {
            os << ' ' << f[fp] + 1;
        }
        os << nl;

    }

    os << "# </faces>" << endl;
}


// ************************************************************************* //
