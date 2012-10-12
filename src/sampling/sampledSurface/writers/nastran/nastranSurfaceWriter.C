/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "nastranSurfaceWriter.H"
#include "IOmanip.H"
#include "Tuple2.H"
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(nastranSurfaceWriter);
    addToRunTimeSelectionTable(surfaceWriter, nastranSurfaceWriter, wordDict);

    // create write methods
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
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nastranSurfaceWriter::formatOS(OFstream& os) const
{
    os.setf(ios_base::scientific);

    // capitalise the E marker
    os.setf(ios_base::uppercase);

    label prec = 0;
    label offset = 7;
    switch (writeFormat_)
    {
        case (wfShort):
        case (wfFree):
        {
            prec = 8 - offset;
            break;
        }
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


void Foam::nastranSurfaceWriter::writeCoord
(
    const point& p,
    label& nPoint,
    label& continuation,
    Ostream& os
) const
{
    // Fixed short/long formats:
    // 1 GRID
    // 2 ID   : point ID
    // 3 CP   : co-ordinate system ID                (blank)
    // 4 X1   : point x cp-ordinate
    // 5 X2   : point x cp-ordinate
    // 6 X3   : point x cp-ordinate
    // 7 CD   : co-ordinate system for displacements (blank)
    // 8 PS   : single point constraints             (blank)
    // 9 SEID : super-element ID

    switch (writeFormat_)
    {
        case wfShort:
        {
            os  << setw(8) << "GRID"
                << setw(8) << ++nPoint
                << "        " 
                << setw(8) << p.x()
                << setw(8) << p.y()
                << setw(8) << p.z()
                << nl;

            break;
        }
        case wfLong:
        {
            os  << setw(8) << "GRID*"
                << setw(16) << ++nPoint
                << "                "
                << setw(16) << p.x()
                << setw(16) << p.y()
                << setw(8) << ++continuation
                << nl
                << setw(8) << continuation
                << setw(16) << p.z()
                << nl;

            break;
        }
        case wfFree:
        {
            os  << "GRID"
                << ',' << ++nPoint
                << ','
                << ',' << p.x()
                << ',' << p.y()
                << ',' << p.z()
                << nl;

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "void Foam::nastranSurfaceWriter::writeCoord"
                "("
                    "Ostream&, "
                    "const point&"
                ") const"
            )   << "Unknown writeFormat enumeration" << abort(FatalError);
        }
    }
}

void Foam::nastranSurfaceWriter::writeFace
(
    const word& faceType,
    const labelList& facePts,
    label& nFace,
    Ostream& os
) const
{
    // Only valid surface elements are CTRIA3 and CQUAD4

    // Fixed short/long formats:
    // 1 CQUAD4
    // 2 EID  : element ID
    // 3 PID  : property element ID; default = EID   (blank)
    // 4 G1   : grid point index
    // 5 G2   : grid point index
    // 6 G3   : grid point index
    // 7 G4   : grid point index
    // 8 onwards - not used

    // For CTRIA3 elements, cols 7 onwards are not used

    switch (writeFormat_)
    {
        case wfShort:
        case wfLong:
        {
            os  << setw(8) << faceType
                << setw(8) << ++nFace
                << "        ";

            forAll(facePts, i)
            {
                os  << setw(8) << facePts[i];
            }

            os  << nl;

            break;
        }
        case wfFree:
        {
            os  << faceType << ','
                << ++nFace << ',';

            forAll(facePts, i)
            {
                os  << ',' << facePts[i];
            }

            os  << nl;

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "void Foam::nastranSurfaceWriter::writeFace"
                "("
                    "const word&"
                    "const labelList&"
                    "label&"
                    "Ostream&, "
                ") const"
            )   << "Unknown writeFormat enumeration" << abort(FatalError);
        }
    }

}


void Foam::nastranSurfaceWriter::writeGeometry
(
    const pointField& points,
    const faceList& faces,
    List<DynamicList<face> >& decomposedFaces,
    Ostream& os
) const
{
    // write points

    os  << "$" << nl
        << "$ Points" << nl
        << "$" << nl;

    label nPoint = 0;
    label continuation = 0;

    forAll(points, pointI)
    {
        writeCoord(points[pointI], nPoint, continuation, os);
    }


    // write faces

    os  << "$" << nl
        << "$ Faces" << nl
        << "$" << nl;

    label nFace = 0;

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
            // decompose poly face into tris
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nastranSurfaceWriter::nastranSurfaceWriter()
:
    surfaceWriter(),
    writeFormat_(wfShort),
    fieldMap_()
{}


Foam::nastranSurfaceWriter::nastranSurfaceWriter(const dictionary& options)
:
    surfaceWriter(),
    writeFormat_(wfShort),
    fieldMap_()
{
    if (options.found("format"))
    {
        writeFormat_ = writeFormatNames_.read(options.lookup("format"));
    }

    List<Tuple2<word, word> > fieldSet(options.lookup("fields"));

    forAll(fieldSet, i)
    {
        fieldMap_.insert(fieldSet[i].first(), fieldSet[i].second());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nastranSurfaceWriter::~nastranSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nastranSurfaceWriter::write
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

    OFstream os(outputDir/surfaceName + ".dat");
    formatOS(os);

    if (verbose)
    {
        Info<< "Writing nastran file to " << os.name() << endl;
    }

    os  << "TITLE=OpeNFOAM " << surfaceName.c_str() << " mesh" << nl
        << "$" << nl
        << "BEGIN BULK" << nl;

    List<DynamicList<face> > decomposedFaces(faces.size());

    writeGeometry(points, faces, decomposedFaces, os);

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    os  << "ENDDATA" << endl;
}


// ************************************************************************* //
