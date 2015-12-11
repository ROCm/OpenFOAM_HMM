/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
#include "Tuple2.H"
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(nastranSurfaceWriter);
    addToRunTimeSelectionTable(surfaceWriter, nastranSurfaceWriter, wordDict);

    // Create write methods
    defineSurfaceWriterWriteFields(nastranSurfaceWriter);

    template<>
    const char* NamedEnum<nastranSurfaceWriter::writeFormat, 3>::names[] =
    {
        "short",
        "long",
        "free"
    };


    const NamedEnum<nastranSurfaceWriter::writeFormat, 3>
        nastranSurfaceWriter::writeFormatNames_;

    template<>
    const char* NamedEnum<nastranSurfaceWriter::dataFormat, 2>::names[] =
    {
        "PLOAD2",
        "PLOAD4"
    };

    const NamedEnum<nastranSurfaceWriter::dataFormat, 2>
        nastranSurfaceWriter::dataFormatNames_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nastranSurfaceWriter::formatOS(OFstream& os) const
{
    os.setf(ios_base::scientific);

    // Capitalise the E marker
    os.setf(ios_base::uppercase);

    label prec = 0;
    label offset = 7;
    switch (writeFormat_)
    {
        case (wfShort):
        {
            prec = 8 - offset;
            break;
        }
        case (wfFree):
        case (wfLong):
        {
            prec = 16 - offset;
            break;
        }
        default:
        {
        }
    }

    os.precision(prec);
}


void Foam::nastranSurfaceWriter::writeKeyword
(
    const word& keyword,
    Ostream& os
) const
{
    os.setf(ios_base::left);

    switch (writeFormat_)
    {
        case wfShort:
        {
            os  << setw(8) << keyword;
            break;
        }
        case wfLong:
        {
            os  << setw(8) << word(keyword + '*');
            break;
        }
        case wfFree:
        {
            os  << keyword;
            break;
        }
    }

    os.unsetf(ios_base::left);
}


void Foam::nastranSurfaceWriter::writeCoord
(
    const point& p,
    const label pointI,
    OFstream& os
) const
{
    // Fixed short/long formats:
    // 1 GRID
    // 2 ID   : point ID - requires starting index of 1
    // 3 CP   : co-ordinate system ID                (blank)
    // 4 X1   : point x cp-ordinate
    // 5 X2   : point x cp-ordinate
    // 6 X3   : point x cp-ordinate
    // 7 CD   : co-ordinate system for displacements (blank)
    // 8 PS   : single point constraints             (blank)
    // 9 SEID : super-element ID


    writeKeyword("GRID", os);

    os  << separator_;

    os.setf(ios_base::right);

    writeValue(pointI + 1, os);
    os  << separator_;
    writeValue("", os);
    os  << separator_;
    writeValue(p.x(), os);
    os  << separator_;
    writeValue(p.y(), os);
    os  << separator_;

    switch (writeFormat_)
    {
        case wfShort:
        {
            os  << setw(8) << p.z()
                << nl;
            os.unsetf(ios_base::right);

            break;
        }
        case wfLong:
        {
            os  << nl;
            os.unsetf(ios_base::right);
            writeKeyword("", os);
            os.setf(ios_base::right);
            writeValue(p.z(), os);
            os  << nl;

            break;
        }
        case wfFree:
        {
            writeValue(p.z(), os);
            os  << nl;

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown writeFormat enumeration" << abort(FatalError);
        }
    }

    os.unsetf(ios_base::right);
}


void Foam::nastranSurfaceWriter::writeFace
(
    const word& faceType,
    const labelList& facePts,
    label& nFace,
    OFstream& os
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

    label PID = 1;

    writeKeyword(faceType, os);

    os  << separator_;

    os.setf(ios_base::right);

    writeValue(nFace++, os);

    os  << separator_;

    writeValue(PID, os);

    switch (writeFormat_)
    {
        case wfShort:
        {
            forAll(facePts, i)
            {
                writeValue(facePts[i] + 1, os);
            }

            break;
        }
        case wfLong:
        {
            forAll(facePts, i)
            {
                writeValue(facePts[i] + 1, os);
                if (i == 1)
                {
                    os  << nl;
                    os.unsetf(ios_base::right);
                    writeKeyword("", os);
                    os.setf(ios_base::right);
                }
            }

            break;
        }
        case wfFree:
        {
            forAll(facePts, i)
            {
                os  << separator_;
                writeValue(facePts[i] + 1, os);
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown writeFormat enumeration" << abort(FatalError);
        }
    }

    os  << nl;
    os.unsetf(ios_base::right);
}


void Foam::nastranSurfaceWriter::writeGeometry
(
    const pointField& points,
    const faceList& faces,
    List<DynamicList<face> >& decomposedFaces,
    OFstream& os
) const
{
    // Write points

    os  << "$" << nl
        << "$ Points" << nl
        << "$" << nl;

    forAll(points, pointI)
    {
        writeCoord(points[pointI], pointI, os);
    }


    // Write faces

    os  << "$" << nl
        << "$ Faces" << nl
        << "$" << nl;

    label nFace = 1;

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        if (f.size() == 3)
        {
            writeFace("CTRIA3", faces[faceI], nFace, os);
            decomposedFaces[faceI].append(faces[faceI]);
        }
        else if (f.size() == 4)
        {
            writeFace("CQUAD4", faces[faceI], nFace, os);
            decomposedFaces[faceI].append(faces[faceI]);
        }
        else
        {
            // Decompose poly face into tris
            label nTri = 0;
            faceList triFaces;
            f.triangles(points, nTri, triFaces);

            forAll(triFaces, triI)
            {
                writeFace("CTRIA3", triFaces[triI], nFace, os);
                decomposedFaces[faceI].append(triFaces[triI]);
            }
        }
    }
}


void Foam::nastranSurfaceWriter::writeFooter(Ostream& os) const
{
    label PID = 1;

    writeKeyword("PSHELL", os);

    os  << separator_;

    writeValue(PID, os);

    for (label i = 0; i < 7; i++)
    {
        // Dummy values
        os  << separator_;
        writeValue(1, os);
    }

    os  << nl;
    writeKeyword("MAT1", os);
    os  << separator_;

    label MID = 1;

    writeValue(MID, os);

    for (label i = 0; i < 7; i++)
    {
        // Dummy values
        os  << separator_;
        writeValue("", os);
    }

    os << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nastranSurfaceWriter::nastranSurfaceWriter()
:
    surfaceWriter(),
    writeFormat_(wfShort),
    fieldMap_(),
    scale_(1.0)
{}


Foam::nastranSurfaceWriter::nastranSurfaceWriter(const dictionary& options)
:
    surfaceWriter(),
    writeFormat_(wfLong),
    fieldMap_(),
    scale_(options.lookupOrDefault("scale", 1.0)),
    separator_("")
{
    if (options.found("format"))
    {
        writeFormat_ = writeFormatNames_.read(options.lookup("format"));
    }

    if (writeFormat_ == wfFree)
    {
        separator_ = ",";
    }

    List<Tuple2<word, word> > fieldSet(options.lookup("fields"));

    forAll(fieldSet, i)
    {
        dataFormat format = dataFormatNames_[fieldSet[i].second()];

        fieldMap_.insert(fieldSet[i].first(), format);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nastranSurfaceWriter::~nastranSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::nastranSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os(outputDir/surfaceName + ".nas");
    formatOS(os);

    if (verbose)
    {
        Info<< "Writing nastran file to " << os.name() << endl;
    }

    os  << "TITLE=OpenFOAM " << surfaceName.c_str() << " mesh" << nl
        << "$" << nl
        << "BEGIN BULK" << nl;

    List<DynamicList<face> > decomposedFaces(faces.size());

    writeGeometry(points, faces, decomposedFaces, os);

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    writeFooter(os);

    os  << "ENDDATA" << endl;

    return os.name();
}


// ************************************************************************* //
