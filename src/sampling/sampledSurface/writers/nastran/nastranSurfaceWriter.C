/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "nastranSurfaceWriter.H"
#include "IOmanip.H"
#include "Pair.H"
#include "HashSet.H"
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(nastranSurfaceWriter);
    addToRunTimeSelectionTable(surfaceWriter, nastranSurfaceWriter, wordDict);

    // Create write methods
    defineSurfaceWriterWriteFields(nastranSurfaceWriter);
}

const Foam::Enum
<
    Foam::nastranSurfaceWriter::loadFormat
>
Foam::nastranSurfaceWriter::loadFormatNames_
{
    { loadFormat::PLOAD2, "PLOAD2" },
    { loadFormat::PLOAD4, "PLOAD4" },
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Ostream& Foam::nastranSurfaceWriter::writeKeyword
(
    Ostream& os,
    const word& keyword
) const
{
    return fileFormats::NASCore::writeKeyword(os, keyword, writeFormat_);
}


void Foam::nastranSurfaceWriter::writeCoord
(
    Ostream& os,
    const point& pt,
    const label pointI
) const
{
    // Fixed short/long formats:
    // 1 GRID
    // 2 ID   : point ID - requires starting index of 1
    // 3 CP   : coordinate system ID                (blank)
    // 4 X1   : point x coordinate
    // 5 X2   : point x coordinate
    // 6 X3   : point x coordinate
    // 7 CD   : coordinate system for displacements (blank)
    // 8 PS   : single point constraints            (blank)
    // 9 SEID : super-element ID

    writeKeyword(os, "GRID")    << separator_;

    os.setf(std::ios_base::right);

    writeValue(os, pointI+1)    << separator_;
    writeValue(os, "")          << separator_;
    writeValue(os, pt.x())      << separator_;
    writeValue(os, pt.y())      << separator_;

    switch (writeFormat_)
    {
        case fieldFormat::SHORT :
        {
            os  << setw(8) << pt.z() << nl;
            os.unsetf(std::ios_base::right);
            break;
        }

        case fieldFormat::LONG :
        {
            os  << nl;
            os.unsetf(std::ios_base::right);
            writeKeyword(os, "");
            os.setf(std::ios_base::right);

            writeValue(os, pt.z())  << nl;
            break;
        }

        case fieldFormat::FREE :
        {
            writeValue(os, pt.z())  << nl;
            break;
        }
    }

    os.unsetf(std::ios_base::right);
}


void Foam::nastranSurfaceWriter::writeFace
(
    Ostream& os,
    const word& faceType,
    const labelUList& facePts,
    const label nFace,
    const label PID
) const
{
    // Only valid surface elements are CTRIA3 and CQUAD4

    // Fixed short/long formats:
    // 1 CQUAD4
    // 2 EID  : element ID
    // 3 PID  : property element ID; default = EID   (blank)
    // 4 G1   : grid point index - requires starting index of 1
    // 5 G2   : grid point index
    // 6 G3   : grid point index
    // 7 G4   : grid point index
    // 8 onwards - not used

    // For CTRIA3 elements, cols 7 onwards are not used

    writeKeyword(os, faceType)   << separator_;

    os.setf(std::ios_base::right);

    writeValue(os, nFace)        << separator_;
    writeValue(os, PID);

    switch (writeFormat_)
    {
        case fieldFormat::SHORT :
        {
            for (const label pointi : facePts)
            {
                writeValue(os, pointi + 1);
            }

            break;
        }

        case fieldFormat::LONG :
        {
            forAll(facePts, i)
            {
                writeValue(os, facePts[i] + 1);
                if (i == 1)
                {
                    os  << nl;
                    os.unsetf(std::ios_base::right);
                    writeKeyword(os, "");
                    os.setf(std::ios_base::right);
                }
            }

            break;
        }

        case fieldFormat::FREE :
        {
            for (const label pointi : facePts)
            {
                os  << separator_;
                writeValue(os, pointi + 1);
            }

            break;
        }
    }

    os  << nl;
    os.unsetf(std::ios_base::right);
}


void Foam::nastranSurfaceWriter::writeGeometry
(
    Ostream& os,
    const meshedSurf& surf,
    List<DynamicList<face>>& decomposedFaces
) const
{
    const pointField& points = surf.points();
    const faceList&    faces = surf.faces();
    const labelList&   zones = surf.zoneIds();

    // Write points

    os  << "$" << nl
        << "$ Points" << nl
        << "$" << nl;

    forAll(points, pointi)
    {
        writeCoord(os, points[pointi], pointi);
    }

    // Write faces
    decomposedFaces.clear();
    decomposedFaces.setSize(faces.size());

    os  << "$" << nl
        << "$ Faces" << nl
        << "$" << nl;

    label nFace = 0; // the element-id
    forAll(faces, facei)
    {
        const face& f = faces[facei];
        // 1-offset for PID
        const label PID = 1 + (facei < zones.size() ? zones[facei] : 0);

        if (f.size() == 3)
        {
            writeFace(os, "CTRIA3", f, ++nFace, PID);
            decomposedFaces[facei].append(f);
        }
        else if (f.size() == 4)
        {
            writeFace(os, "CQUAD4", f, ++nFace, PID);
            decomposedFaces[facei].append(f);
        }
        else
        {
            // Decompose poly face into tris
            label nTri = 0;
            faceList triFaces;
            f.triangles(points, nTri, triFaces);

            forAll(triFaces, trii)
            {
                writeFace(os, "CTRIA3", triFaces[trii], ++nFace, PID);
                decomposedFaces[facei].append(triFaces[trii]);
            }
        }
    }
}


Foam::Ostream& Foam::nastranSurfaceWriter::writeFooter
(
    Ostream& os,
    const meshedSurf& surf
) const
{
    // zone id have been used for the PID. Find unique values.

    labelList pidsUsed = labelHashSet(surf.zoneIds()).sortedToc();
    if (pidsUsed.empty())
    {
        pidsUsed.setSize(1, Zero); // fallback
    }

    for (auto pid : pidsUsed)
    {
        writeKeyword(os, "PSHELL")  << separator_;
        writeValue(os, pid+1);  // 1-offset for PID

        for (label i = 0; i < 7; ++i)
        {
            // Dummy values
            os  << separator_;
            writeValue(os, 1);
        }
        os  << nl;
    }


    // use single material ID

    const label MID = 1;

    writeKeyword(os, "MAT1")    << separator_;
    writeValue(os, MID);

    for (label i = 0; i < 7; ++i)
    {
        // Dummy values
        os  << separator_;
        writeValue(os, "");
    }
    os << nl;

    return os;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nastranSurfaceWriter::nastranSurfaceWriter()
:
    surfaceWriter(),
    writeFormat_(fieldFormat::SHORT),
    fieldMap_(),
    scale_(1.0),
    separator_()
{}


Foam::nastranSurfaceWriter::nastranSurfaceWriter(const dictionary& options)
:
    surfaceWriter(),
    writeFormat_
    (
        fileFormats::NASCore::fieldFormatNames.lookupOrDefault
        (
            "format",
            options,
            fieldFormat::LONG
        )
    ),
    fieldMap_(),
    scale_(options.lookupOrDefault<scalar>("scale", 1.0)),
    separator_()
{
    if (writeFormat_ == fieldFormat::FREE)
    {
        separator_ = ",";
    }

    List<Pair<word>> fieldPairs(options.lookup("fields"));

    for (const Pair<word>& item : fieldPairs)
    {
        // (field name => load format)
        fieldMap_.insert
        (
            item.first(),
            loadFormatNames_[item.second()]
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::nastranSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os(outputDir/surfaceName + ".nas");
    fileFormats::NASCore::setPrecision(os, writeFormat_);

    if (verbose)
    {
        Info<< "Writing nastran file to " << os.name() << endl;
    }

    os  << "TITLE=OpenFOAM " << surfaceName.c_str()
        << " mesh" << nl
        << "$" << nl
        << "BEGIN BULK" << nl;

    List<DynamicList<face>> decomposedFaces;
    writeGeometry(os, surf, decomposedFaces);

    writeFooter(os, surf)
        << "ENDDATA" << nl;

    return os.name();
}


// ************************************************************************* //
