/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "x3dSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "MeshedSurfaceProxy.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(x3dWriter);
    addToRunTimeSelectionTable(surfaceWriter, x3dWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, x3dWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- A (0-1) range for colouring
template<class Type>
static inline scalar rangex(const scalarMinMax& range, const Type& val)
{
    scalar x = Foam::mag(val);

    return (x - range.min()) / (range.max() - range.min());
}


//- A (0-1) range for colouring
template<>
inline scalar rangex(const scalarMinMax& range, const scalar& val)
{
    scalar x = val;
    return (x - range.min()) / (range.max() - range.min());
}


//- A (0-1) range for colouring
template<>
inline scalar rangex(const scalarMinMax& range, const label& val)
{
    scalar x = val;
    return (x - range.min()) / (range.max() - range.min());
}


static inline void printColour(Ostream& os, const vector& rgb)
{
    os  << rgb[0] << ' ' << rgb[1] << ' ' << rgb[2] << ',' << nl;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::x3dWriter::x3dWriter()
:
    surfaceWriter(),
    range_(),
    colourTablePtr_(nullptr)
{}


Foam::surfaceWriters::x3dWriter::x3dWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    range_(),
    colourTablePtr_(nullptr)
{
    verbose_ = true;

    options.readIfPresent("range", range_);

    word tableName;
    if (options.readIfPresent("colourMap", tableName))
    {
        colourTablePtr_ = colourTable::ptr(tableName);
        if (!colourTablePtr_)
        {
            WarningInFunction
                << "No colourMap " << tableName << " using default" << nl;
        }
    }

    if (!colourTablePtr_)
    {
        tableName = colourTable::predefinedNames[colourTable::COOL_WARM];
        colourTablePtr_ = colourTable::ptr(colourTable::COOL_WARM);
    }

    if (verbose_)
    {
        Info<< "X3D with colourMap '" << tableName << "' and range ";

        if (range_.valid())
        {
            Info<< range_;
        }
        else
        {
            Info<< "auto";
        }
        Info<< nl;
    }
}


Foam::surfaceWriters::x3dWriter::x3dWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    x3dWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::x3dWriter::x3dWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    x3dWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::x3dWriter::write()
{
    checkOpen();

    // Geometry:  rootdir/<TIME>/surfaceName.x3d

    fileName outputFile = outputPath_;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile = outputPath_.path() / timeName() / outputPath_.name();
    }
    outputFile.ext("x3d");

    if (verbose_)
    {
        Info<< "Writing geometry to " << outputFile << endl;
    }

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        MeshedSurfaceProxy<face>
        (
            surf.points(),
            surf.faces()
        ).write(outputFile, "x3d");
    }

    wroteGeom_ = true;
    return outputFile;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::surfaceWriters::x3dWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    if (!colourTablePtr_)
    {
        // Write geometry only if there are no colours to use
        WarningInFunction
            << "No output colours set" << endl;

        return this->write();
    }

    checkOpen();

    // Field:  rootdir/<TIME>/<field>_surfaceName.x3d

    fileName outputFile = outputPath_.path();
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile /= timeName();
    }

    // Append <field>_surfaceName.usr
    outputFile /= fieldName + '_' + outputPath_.name();
    outputFile.ext("x3d");

    if (verbose_)
    {
        Info<< "Writing field " << fieldName << " to " << outputFile << endl;
    }

    const meshedSurf& surf = surface();

    // geometry merge() implicit
    tmp<Field<Type>> tfield = mergeField(localValues);

    if (Pstream::master() || !parallel_)
    {
        const auto& values = tfield();

        scalarMinMax range(range_);

        if (!range.valid())
        {
            range = minMaxMag(values);
        }

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        OFstream os(outputFile);

        writeHeader(os);
        beginGroup(os);
        writeAppearance(os);

        // For point field: "colorPerVetex=true"
        os  << "  <IndexedFaceSet"
            << " colorPerVertex='" << Switch(this->isPointData()) << "'"
            << " coordIndex='" << nl;

        for (const auto& f : surf.faces())
        {
            for (const label vrti : f)
            {
                os << vrti << ' ';
            }
            os << "-1\n";
        }
        os  << "'";

        // Colour indices for face fields
        if (!this->isPointData())
        {
            const label nFaces = surf.faces().size();

            os  << " colorIndex='";

            for (label i=0; i < nFaces; ++i)
            {
                os << i << ' ';
            }
            os  << "'";
        }

        os  << " >\n";  // IndexedFaceSet

        writePoints(os, surf.points());

        os << "<Color color='" << nl;

        // writeColours(os, values, range, colorBar);

        for (const Type& val : values)
        {
            const scalar x = rangex(range, val);
            vector rgb = colourTablePtr_->value(x);
            printColour(os, rgb);
        }

        os << "' />" << nl;  // Color

        os  <<
            "   </IndexedFaceSet>\n";

        endGroup(os);
        writeFooter(os);
    }

    wroteGeom_ = true;
    return outputFile;
}


// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::x3dWriter);


// ************************************************************************* //
