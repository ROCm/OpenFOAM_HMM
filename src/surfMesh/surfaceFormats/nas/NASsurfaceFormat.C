/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "NASsurfaceFormat.H"
#include "IFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::NASsurfaceFormat<Face>::NASsurfaceFormat()
:
    ParentType()
{}


template<class Face>
Foam::fileFormats::NASsurfaceFormat<Face>::NASsurfaceFormat
(
    const fileName& fName
)
:
    ParentType()
{
    ThisType::read(fName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::NASsurfaceFormat<Face>::read
(
    const fileName& fName
)
{
    ParentType::clear();

    // triangulation required?
    bool mustTriangulate = false;
    {
        Face f;
        if (f.max_size() == 3)
        {
            mustTriangulate = true;
        }
    }

    IFstream is(fName);
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::NASsurfaceFormat::read(const fileName&)"
        )
            << "Cannot read file " << fName
            << exit(FatalError);
    }

    DynamicList<point>  pointLst;
    // Nastran index of points
    DynamicList<label>  pointId;
    DynamicList<Face>   faceLst;
    DynamicList<label>  regionLst;
    HashTable<label>    groupToPatch;

    // From face groupId to patchId
    Map<label> groupIdToPatchId;
    label nPatches = 0;

    // Name for face group
    Map<word> groupIdToName;

    // Ansa tags. Denoted by $ANSA_NAME. These will appear just before the
    // first use of a type. We read them and store the pshell types which
    // are used to name the patches.
    label ansaId = -1;
    word ansaType;
    word ansaName;

    // leave faces that didn't have a group in 0
    //    label groupID = 0;
    //    label maxGroupID = -1;

    // A single warning per unrecognized command
    HashSet<word> unhandledCmd;

    while (is.good())
    {
        string line;
        is.getLine(line);

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
                ansaId = readLabel
                (
                    IStringStream(line.substr(sem0+1, sem1-sem0-1))()
                );
                ansaType = line.substr(sem1+1, sem2-sem1-1);

                string rawName;
                is.getLine(rawName);
                if (rawName[rawName.size()-1] == '\r')
                {
                    rawName = rawName.substr(1, rawName.size()-2);
                }
                else
                {
                    rawName = rawName.substr(1, rawName.size()-1);
                }

                string::stripInvalid<word>(rawName);
                ansaName = rawName;

                // Info<< "ANSA tag for NastranID:" << ansaID
                //     << " of type " << ansaType
                //     << " name " << ansaName << endl;
            }
        }


        // Hypermesh extension
        // $HMNAME COMP                   1"partName"
        if
        (
            line.substr(0, 12) == "$HMNAME COMP"
         && line.find ('"') != string::npos
        )
        {
            label groupId = readLabel
            (
                IStringStream(line.substr(16, 16))()
            );

            IStringStream lineStream(line.substr(32));

            string rawName;
            lineStream >> rawName;
            string::stripInvalid<word>(rawName);

            word groupName(rawName);
            groupIdToName.insert(groupId, groupName);

            Info<< "group " << groupId << " => " << groupName << endl;
        }


        // Skip empty or comment
        if (line.size() == 0 || line[0] == '$')
        {
            continue;
        }


        // Check if character 72 is continuation
        if (line.size() > 72 && line[72] == '+')
        {
            line = line.substr(0, 72);

            while (true)
            {
                string buf;
                is.getLine(buf);

                if (buf.size() > 72 && buf[72] == '+')
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
            triFace fTri;

            label groupId = readLabel(IStringStream(line.substr(16,8))());
            fTri[0] = readLabel(IStringStream(line.substr(24,8))());
            fTri[1] = readLabel(IStringStream(line.substr(32,8))());
            fTri[2] = readLabel(IStringStream(line.substr(40,8))());

            // Convert groupID into patchID
            Map<label>::const_iterator iter = groupIdToPatchId.find(groupId);

            label patchI;
            if (iter == groupIdToPatchId.end())
            {
                patchI = nPatches++;
                groupIdToPatchId.insert(groupId, patchI);
                Info<< "patch " << patchI << " => group " << groupId << endl;
            }
            else
            {
                patchI = iter();
            }

            faceLst.append(fTri);
            regionLst.append(patchI);
        }
        else if (cmd == "CQUAD4")
        {
            face fQuad(4);
            UList<label>& f = static_cast<UList<label>&>(fQuad);

            label groupId = readLabel(IStringStream(line.substr(16,8))());
            fQuad[0] = readLabel(IStringStream(line.substr(24,8))());
            fQuad[1] = readLabel(IStringStream(line.substr(32,8))());
            fQuad[2] = readLabel(IStringStream(line.substr(40,8))());
            fQuad[3] = readLabel(IStringStream(line.substr(48,8))());

            // Convert group into patch
            Map<label>::const_iterator iter = groupIdToPatchId.find(groupId);

            label patchI;
            if (iter == groupIdToPatchId.end())
            {
                patchI = nPatches++;
                groupIdToPatchId.insert(groupId, patchI);
                Info<< "patch " << patchI << " => group " << groupId << endl;
            }
            else
            {
                patchI = iter();
            }

            if (mustTriangulate)
            {
                triFace fTri;

                // simple face triangulation about f[0].
                // cannot use face::triangulation since points are incomplete
                fTri[0] = f[0];
                for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
                {
                    label fp2 = (fp1 + 1) % f.size();

                    fTri[1] = f[fp1];
                    fTri[2] = f[fp2];

                    faceLst.append(fTri);
                    regionLst.append(patchI);
                }
            }
            else
            {
                faceLst.append(Face(f));
                regionLst.append(patchI);
            }
        }
        else if (cmd == "PSHELL")
        {
            // Read shell type since group gives patchnames.
            label groupId = readLabel(IStringStream(line.substr(8,8))());

            if (groupId == ansaId && ansaType == "PSHELL")
            {
                groupIdToName.insert(groupId, ansaName);
                Info<< "group " << groupId << " => " << ansaName << endl;
            }
        }
        else if (cmd == "GRID")
        {
            label index = readLabel(IStringStream(line.substr(8,8))());
            scalar x = parseNASCoord(line.substr(24, 8));
            scalar y = parseNASCoord(line.substr(32, 8));
            scalar z = parseNASCoord(line.substr(40, 8));

            pointId.append(index);
            pointLst.append(point(x, y, z));
        }
        else if (cmd == "GRID*")
        {
            // Long format is on two lines with '*' continuation symbol
            // on start of second line.
            // Typical line (spaces compacted)
            // GRID*      126   0 -5.55999875E+02 -5.68730474E+02
            // *         2.14897901E+02

            label index = readLabel(IStringStream(line.substr(8,16))());
            scalar x = parseNASCoord(line.substr(40, 16));
            scalar y = parseNASCoord(line.substr(56, 16));

            is.getLine(line);
            if (line[0] != '*')
            {
                FatalErrorIn
                (
                    "fileFormats::NASsurfaceFormat::read(const fileName&)"
                )
                    << "Expected continuation symbol '*' when reading GRID*"
                    << " (double precision coordinate) output" << nl
                    << "Read:" << line << nl
                    << "File:" << is.name()
                    << " line:" << is.lineNumber()
                    << exit(FatalError);
            }
            scalar z = parseNASCoord(line.substr(8, 16));

            pointId.append(index);
            pointLst.append(point(x, y, z));
        }
        else if (unhandledCmd.insert(cmd))
        {
            Info<< "Unhandled Nastran command " << line << nl
                << "File:" << is.name()
                << " line:" << is.lineNumber()
                << endl;
        }
    }

    Info<< "Read faces:" << faceLst.size()
        << " points:" << pointLst.size()
        << endl;


    // transfer to normal lists
    ParentType::points().transfer(pointLst);
    ParentType::regions().transfer(regionLst);

    pointId.shrink();
    faceLst.shrink();

    {
        // Build inverse mapping (index to point)
        Map<label> nasToFoamPoint(2*pointId.size());
        forAll(pointId, i)
        {
            nasToFoamPoint.insert(pointId[i], i);
        }
        pointId.clearStorage();


        // Relabel faces
        forAll(faceLst, i)
        {
            Face& f = faceLst[i];
            forAll(f, fp)
            {
                f[fp] = nasToFoamPoint[f[fp]];
            }
        }
    }


    // convert Nastran groupId => name to patchId => name
    Map<word> regionNames;
    forAllConstIter(Map<word>, groupIdToName, iter)
    {
        Map<label>::const_iterator iter2 = groupIdToPatchId.find(iter.key());

        if (iter2 != groupIdToPatchId.end())
        {
            regionNames.insert(iter2.key(), iter());
        }
    }

    // transfer to normal lists
    ParentType::faces().transfer(faceLst);

    ParentType::setPatches(regionNames);
    return true;
}

// ************************************************************************* //
