/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
#include "ensightOutput.H"
#include "ensightOutputSurface.H"
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
    const scalar timeValue
)
{
    os
        << "time set:               " << ts << nl
        << "number of steps:        " << 1 << nl;

    // Single value - starts at index 0
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
    label pos_;

    os
        << "time set:               " << ts << nl
        << "number of steps:        " << values.size() << nl;

    // Assume contiguous numbering - starts at index 0
    os  << "filename start number:  0" << nl
        << "filename increment:     1" << nl;


    os  << "time values:" << nl;
    pos_ = 0;
    for (const scalar& val : values)
    {
        if (pos_ == 6)
        {
            os  << nl;
            pos_ = 0;
        }
        ++pos_;

        os  << ' ' << setf(ios_base::right) << setw(12) << val;
    }
    os  << nl << nl;
}


void Foam::surfaceWriters::ensightWriter::printTimeset
(
    OSstream& os,
    const label ts,
    const UList<scalar>& values,
    const bitSet& indices
)
{
    label pos_;

    // Check if continuous numbering can be used
    if
    (
        values.empty()
     || (indices.size() == values.size() && indices.all())
    )
    {
        // Can simply emit as 0-based with increment
        printTimeset(os, ts, values);
        return;
    }


    // Generate time set
    os
        << "time set:               " << ts << nl
        << "number of steps:        " << indices.count() << nl;


    os  << "filename numbers:" << nl;
    pos_ = 0;
    for (const label& idx : indices)
    {
        if (pos_ == 6)
        {
            os  << nl;
            pos_ = 0;
        }
        ++pos_;

        os  << ' ' << setf(ios_base::right) << setw(8) << idx;
    }
    os  << nl;


    os  << "time values:" << nl;
    pos_ = 0;
    for (const label& idx : indices)
    {
        if (pos_ == 6)
        {
            os  << nl;
            pos_ = 0;
        }
        ++pos_;

        os  << ' ' << setf(ios_base::right) << setw(12) << values[idx];
    }
    os  << nl << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::ensightWriter::ensightWriter()
:
    surfaceWriter(),
    writeFormat_(IOstream::ASCII),
    collateTimes_(true),
    caching_("fieldsDict")  // Historic name
{}


Foam::surfaceWriters::ensightWriter::ensightWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    writeFormat_
    (
        IOstreamOption::formatEnum("format", options, IOstream::ASCII)
    ),
    collateTimes_(options.getOrDefault("collateTimes", true)),
    caching_("fieldsDict")  // Historic name
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

void Foam::surfaceWriters::ensightWriter::close()
{
    caching_.clear();
    surfaceWriter::close();
}


// Note that ensight does supports geometry in a separate file,
// but setting this true leaves geometry files in the wrong places
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
