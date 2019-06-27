/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "ensightSurfaceWriter.H"
#include "IOmanip.H"
#include "Fstream.H"
#include "OSspecific.H"
#include "ensightCase.H"
#include "ensightPartFaces.H"
#include "ensightOutput.H"
#include "ensightPTraits.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(ensightWriter);
    addToRunTimeSelectionTable(surfaceWriter, ensightWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, ensightWriter, wordDict);
}
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::surfaceWriters::ensightWriter::printTimeset
(
    OSstream& os,
    const label ts,
    const scalar& timeValue
)
{
    os
        << "time set:               " << ts << nl
        << "number of steps:        " << 1 << nl;

    // Assume to be contiguous numbering
    os  << "filename start number:  0" << nl
        << "filename increment:     1" << nl
        << "time values:" << nl;

    os  << "    " << timeValue
        << nl << nl;
}


void Foam::surfaceWriters::ensightWriter::printTimeset
(
    OSstream& os,
    const label ts,
    const UList<scalar>& values
)
{
    label count = values.size();
    os
        << "time set:               " << ts << nl
        << "number of steps:        " << count << nl;

    // Assume to be contiguous numbering
    os  << "filename start number:  0" << nl
        << "filename increment:     1" << nl
        << "time values:" << nl;

    count = 0;
    for (const scalar& t : values)
    {
        os << ' ' << setw(12) << t;

        if (++count % 6 == 0)
        {
            os << nl;
        }
    }
    os  << nl << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::ensightWriter::ensightWriter()
:
    surfaceWriter(),
    writeFormat_(IOstream::ASCII),
    collateTimes_(true)
{}


Foam::surfaceWriters::ensightWriter::ensightWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    writeFormat_
    (
        IOstreamOption::formatNames.lookupOrDefault
        (
            "format",
            options,
            IOstreamOption::ASCII,
            true  // Failsafe behaviour
        )
    ),
    collateTimes_(options.lookupOrDefault("collateTimes", true))
{}


Foam::surfaceWriters::ensightWriter::ensightWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    ensightWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::ensightWriter::ensightWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    ensightWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Note that ensight does supports geometry in a separate file,
// but setting this true leaves mesh files in the wrong places
// (when there are fields).
//
// Make this false to let the field writers take back control
bool Foam::surfaceWriters::ensightWriter::separateGeometry() const
{
    return false;
}


Foam::fileName Foam::surfaceWriters::ensightWriter::write()
{
    // if (collateTimes_)
    // {
    //     return writeCollated();
    // }
    // else
    // {
    //     return writeUncollated();
    // }

    return writeUncollated();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightSurfaceWriterCollated.C"
#include "ensightSurfaceWriterUncollated.C"


// Field writing implementations

template<class Type>
Foam::fileName Foam::surfaceWriters::ensightWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    if (collateTimes_)
    {
        return writeCollated(fieldName, localValues);
    }
    else
    {
        return writeUncollated(fieldName, localValues);
    }
}


defineSurfaceWriterWriteFields(Foam::surfaceWriters::ensightWriter);


// ************************************************************************* //
