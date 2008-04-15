/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description
    Nastran surface reader. Does Ansa $ANSA_NAME extension to get name
    of patch. Handles Ansa coordinates like:

        GRID          28        10.20269-.030265-2.358-8

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "IFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Do weird things to extract number
static scalar parseNASCoord(const string& s)
{
    size_t expSign = s.find_last_of("+-");

    if (expSign != string::npos && expSign > 0 && !isspace(s[expSign-1]))
    {
        scalar mantissa = readScalar(IStringStream(s.substr(0, expSign))());
        scalar exp = readScalar(IStringStream(s.substr(expSign+1))());

        if (s[expSign] == '-')
        {
            exp = -exp;
        }
        return mantissa*pow(10, exp);
    }
    else
    {
        return readScalar(IStringStream(s)());
    }
}


bool triSurface::readNAS(const fileName& OBJfileName)
{
    IFstream OBJfile(OBJfileName);

    if (!OBJfile.good())
    {
        FatalErrorIn("triSurface::readNAS(const fileName&)")
            << "Cannot read file " << OBJfileName
            << exit(FatalError);
    }

    // coordinates of point
    DynamicList<point> points;
    // Nastran index of point
    DynamicList<label> indices;
    // Faces in terms of Nastran point indices
    DynamicList<labelledTri> faces;
    // From face group to patch
    Map<label> groupToPatch;
    label nPatches = 0;
    // Name for face group
    Map<word> groupToName;

    // Ansa tags. Denoted by $ANSA_NAME. These will appear just before the
    // first use of a type. We read them and store the pshell types which
    // are used to name the patches.
    label ansaID = -1;
    word ansaType;
    string ansaName;

    // Done warnings per unrecognized command
    HashSet<word> unhandledCmd;

    while (OBJfile.good())
    {
        string line;
        OBJfile.getLine(line);

        // Ansa extension
        if (line.substr(0, 10) == "$ANSA_NAME")
        {
            string::size_type sem0 = line.find (';', 0);
            string::size_type sem1 = line.find (';', sem0+1);
            string::size_type sem2 = line.find (';', sem1+1);

            if
            (
                sem0 != string::npos
             && sem1 != string::npos
             && sem2 != string::npos
            )
            {
                ansaID = readLabel
                (
                    IStringStream(line.substr(sem0+1, sem1-sem0-1))()
                );
                ansaType = line.substr(sem1+1, sem2-sem1-1);

                string nameString;
                OBJfile.getLine(ansaName);
                if (ansaName[ansaName.size()-1] == '\r')
                {
                    ansaName = ansaName.substr(1, ansaName.size()-2);
                }
                else
                {
                    ansaName = ansaName.substr(1, ansaName.size()-1);
                }
                //Pout<< "ANSA tag for NastranID:" << ansaID
                //    << " of type " << ansaType
                //    << " name " << ansaName << endl;
            }
        }


        if (line.size() == 0 || line[0] == '$')
        {
            // Skip empty or comment
            continue;
        }

        // Check if character 72 is continuation
        if (line.size() > 72 && line[72] == '+')
        {
            line = line.substr(0, 72);

            while (true)
            {
                string buf;
                OBJfile.getLine(buf);

                if (buf.size() > 72 && buf[72]=='+')
                {
                    line += buf.substr(8, 64);
                }
                else
                {
                    line += buf.substr(8, buf.size()-8);
                    break;
                }
            }
        }

        // Read first word
        IStringStream lineStream(line);
        word cmd;
        lineStream >> cmd;

        if (cmd == "CTRIA3")
        {
            //label index, group, a, b, c;
            //lineStream >> index >> group >> a >> b >> c;
            label group = readLabel(IStringStream(line.substr(16,8))());
            label a = readLabel(IStringStream(line.substr(24,8))());
            label b = readLabel(IStringStream(line.substr(32,8))());
            label c = readLabel(IStringStream(line.substr(40,8))());


            // Convert group into patch
            Map<label>::const_iterator iter = groupToPatch.find(group);

            label patchI;
            if (iter == groupToPatch.end())
            {
                patchI = nPatches++;

                Pout<< "Allocating Foam patch " << patchI
                    << " for group " << group << endl;

                groupToPatch.insert(group, patchI);
            }
            else
            {
                patchI = iter();
            }

            faces.append(labelledTri(a, b, c, patchI));
        }
        else if (cmd == "CQUAD4")
        {
            //label index, group, a, b, c, d;
            //lineStream >> index >> group >> a >> b >> c >> d;
            label group = readLabel(IStringStream(line.substr(16,8))());
            label a = readLabel(IStringStream(line.substr(24,8))());
            label b = readLabel(IStringStream(line.substr(32,8))());
            label c = readLabel(IStringStream(line.substr(40,8))());
            label d = readLabel(IStringStream(line.substr(48,8))());

            // Convert group into patch
            Map<label>::const_iterator iter = groupToPatch.find(group);

            label patchI;
            if (iter == groupToPatch.end())
            {
                patchI = nPatches++;

                Pout<< "Allocating Foam patch " << patchI
                    << " for group " << group << endl;

                groupToPatch.insert(group, patchI);
            }
            else
            {
                patchI = iter();
            }

            faces.append(labelledTri(a, b, c, patchI));
            faces.append(labelledTri(c, d, a, patchI));
        }
        else if (cmd == "PSHELL")
        {
            // Read shell type since gives patchnames.
            //label group;
            //lineStream >> group;
            label group = readLabel(IStringStream(line.substr(8,8))());

            if (group == ansaID && ansaType == "PSHELL")
            {
                Pout<< "Found name " << ansaName << " for group "
                    << group << endl;
                groupToName.insert(group, string::validate<word>(ansaName));
            }
        }
        else if (cmd == "GRID")
        {
            //label index;
            //lineStream >> index;
            label index = readLabel(IStringStream(line.substr(8,8))());
            indices.append(index);

            scalar x = parseNASCoord(line.substr(24, 8));
            scalar y = parseNASCoord(line.substr(32, 8));
            scalar z = parseNASCoord(line.substr(40, 8));
            points.append(point(x, y, z));
        }
        else if (cmd == "GRID*")
        {
            // Assume on two lines with '*' continuation symbol on start of
            // second line. (comes out of Tgrid. Typical line (spaces truncated)
            // GRID*      126   0 -5.55999875E+02 -5.68730474E+02
            // *         2.14897901E+02
            string line2;
            OBJfile.getLine(line2);
            if (line2[0] != '*')
            {
                FatalErrorIn("triSurface::readNAS(const fileName&)")
                    << "Expected continuation symbol '*' when reading GRID*"
                    << " (double precision coordinate) output by Tgrid" << nl
                    << "Read:" << line2 << nl
                    << "File:" << OBJfile.name()
                    << " line:" << OBJfile.lineNumber()
                    << exit(FatalError);
            }
            IStringStream lineStream(line.substr(10) + line2.substr(1));

            label index;
            lineStream >> index;
            indices.append(index);

            readScalar(lineStream); // What is this field?
            scalar x = readScalar(lineStream);
            scalar y = readScalar(lineStream);
            scalar z = readScalar(lineStream);
            points.append(point(x, y, z));
        }
        else if (unhandledCmd.insert(cmd))
        {
            Info<< "Unhandled Nastran command " << line << nl
                << "File:" << OBJfile.name()
                << " line:" << OBJfile.lineNumber()
                << endl;
        }
    }

    points.shrink();
    indices.shrink();
    faces.shrink();


    Pout<< "Read triangles:" << faces.size() << " points:" << points.size()
        << endl;

    {
        // Build inverse mapping (index to point)
        Map<label> indexToPoint(2*indices.size());
        forAll(indices, i)
        {
            indexToPoint.insert(indices[i], i);
        }

        // Relabel triangles
        forAll(faces, i)
        {
            labelledTri& f = faces[i];

            f[0] = indexToPoint[f[0]];
            f[1] = indexToPoint[f[1]];
            f[2] = indexToPoint[f[2]];
        }
    }


    // Convert groupToPatch to patchList.
    geometricSurfacePatchList patches(nPatches);

    forAllConstIter(Map<word>, groupToName, iter)
    {
        label patchI = groupToPatch[iter.key()];

        patches[patchI] = geometricSurfacePatch
        (
            "empty",
            iter(),
            patchI
        );
    }

    Pout<< "patches:" << patches << endl;


    // Transfer DynamicLists to straight ones.
    pointField allPoints;
    allPoints.transfer(points);
    points.clear();

    // Create triSurface
    *this = triSurface(faces, patches, allPoints);

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
