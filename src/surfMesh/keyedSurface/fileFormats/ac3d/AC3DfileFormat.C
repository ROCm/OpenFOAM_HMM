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

#include "AC3DfileFormat.H"
#include "triFace.H"
#include "clock.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "tensor.H"
#include "triFace.H"
#include "primitivePatch.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace fileFormats
{

addNamedToRunTimeSelectionTable
(
    keyedSurface,
    AC3DfileFormat,
    fileExtension,
    ac
);

addNamedToMemberFunctionSelectionTable
(
    keyedSurface,
    AC3DfileFormat,
    write,
    fileExtension,
    ac
);

addNamedToMemberFunctionSelectionTable
(
    meshedSurface,
    AC3DfileFormat,
    write,
    fileExtension,
    ac
);

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Type Foam::fileFormats::AC3DfileFormat::parse(const string& str)
{
    IStringStream ss(str);

    Type t;
    ss >> t;
    return t;
}


bool Foam::fileFormats::AC3DfileFormat::readCmd
(
    IFstream& is,
    string& cmd,
    string& args
)
{
    if (is.good())
    {
        string line;
        is.getLine(line);

        string::size_type space = line.find(' ');

        if (space != string::npos)
        {
            cmd  = line.substr(0, space);
            args = line.substr(space+1);

            return true;
        }
    }
    return false;
}


// Read up to line starting with cmd. Sets args to rest of line.
// Returns true if found, false if stream is not good anymore.
bool Foam::fileFormats::AC3DfileFormat::cueTo
(
    IFstream& is,
    const string& cmd,
    string& args
)
{
    while (is.good())
    {
        string line;
        is.getLine(line);

        string::size_type space = line.find(' ');

        if (space != string::npos)
        {
            if (line.substr(0, space) == cmd)
            {
                args = line.substr(space+1);

                return true;
            }
        }
    }
    return false;
}


// Similar to cueTo(), but throws error if cmd not found
Foam::string Foam::fileFormats::AC3DfileFormat::cueToOrDie
(
    IFstream& is,
    const string& cmd,
    const string& errorMsg
)
{
    string args;
    if (!cueTo(is, cmd, args))
    {
        FatalErrorIn
        (
            "fileFormats::AC3DfileFormat::AC3DfileFormat"
            "(const fileName&)"
        )
            << "Cannot find command " << cmd
            << " " << errorMsg
            << exit(FatalError);
    }

    return args;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::AC3DfileFormat::AC3DfileFormat()
:
    Foam::keyedSurface()
{}


Foam::fileFormats::AC3DfileFormat::AC3DfileFormat
(
    const fileName& fName,
    const bool triangulate
)
:
    Foam::keyedSurface()
{
    IFstream is(fName);

    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::AC3DfileFormat::AC3DfileFormat"
            "(const fileName&)"
        )
            << "Cannot read file " << fName
            << exit(FatalError);
    }

    string line, cmd, args;

    is.getLine(line);

    string version = line.substr(4);

    if (version != "b")
    {
        WarningIn
        (
            "fileFormats::AC3DfileFormat::AC3DfileFormat"
            "(const fileName&)"
        )
            << "When reading AC3D file " << fName
            << " read header " << line << " with version "
            << version << endl
            << "Only tested reading with version 'b'."
            << " This might give problems" << endl;
    }


    if (!cueTo(is, "OBJECT", args) || (args != "world"))
    {
        FatalErrorIn
        (
            "fileFormats::AC3DfileFormat::AC3DfileFormat"
            "(const fileName&)"
        )
            << "Cannot find \"OBJECT world\" in file " << fName
            << exit(FatalError);
    }

    // # of kids is the # of patches
    args = cueToOrDie(is, "kids");
    label nPatches = parse<int>(args);

    // Start of vertices for object/patch
    label patchVertOffset = 0;

    DynamicList<point>    pointLst;
    DynamicList<FaceType> faceLst;
    DynamicList<label>    regionLst;

    // patchId => patchName
    Map<word> regionNames;

    for (label patchI = 0; patchI < nPatches; ++patchI)
    {
        word patchName = word("patch") + Foam::name(patchI);

        args = cueToOrDie(is, "OBJECT", "while reading " + patchName);

        // number of vertices for this patch
        label  nPatchPoints = 0;
        vector location(pTraits<vector>::zero);
        // tensor rotation(I);

        // Read all info for current patch
        while (is.good())
        {
            // Read line and get first word. If end of file break since
            // patch should always end with 'kids' command ?not sure.
            if (!readCmd(is, cmd, args))
            {
                FatalErrorIn
                (
                    "fileFormats::AC3DfileFormat::AC3DfileFormat"
                    "(const fileName&)"
                )
                    << "Did not read up to \"kids 0\" while reading patch "
                    << patchI << " from file " << fName
                    << exit(FatalError);
            }

            if (cmd == "name")
            {
                // name %s
                string str = parse<string>(args);
                string::stripInvalid<word>(str);

                patchName = str;
            }
            else if (cmd == "rot")
            {
                // rot  %f %f %f  %f %f %f  %f %f %f

                // IStringStream lineStream(args);
                //
                // lineStream
                //     >> rotation.xx() >> rotation.xy() >> rotation.xz()
                //     >> rotation.yx() >> rotation.yy() >> rotation.yz()
                //     >> rotation.zx() >> rotation.zy() >> rotation.zz();

                WarningIn
                (
                    "fileFormats::AC3DfileFormat::AC3DfileFormat"
                    "(const fileName&)"
                )
                    << "rot (rotation tensor) command not implemented"
                    << "Line:" << cmd << ' ' << args << endl
                    << "while reading patch " << patchI << endl;
            }
            else if (cmd == "loc")
            {
                // loc  %f %f %f
                IStringStream lineStream(args);

                lineStream
                    >> location.x()
                    >> location.y()
                    >> location.z();
            }
            else if (cmd == "numvert")
            {
                // numvert  %d
                nPatchPoints = parse<int>(args);

                for (label vertI = 0; vertI < nPatchPoints; ++vertI)
                {
                    is.getLine(line);
                    IStringStream lineStream(line);

                    point pt;
                    lineStream
                        >> pt.x() >> pt.y() >> pt.z();

                    // Offset with current translation vector
                    pointLst.append(location + pt);
                }
            }
            else if (cmd == "numsurf")
            {
                label nFaces = parse<int>(args);

                for (label faceI = 0; faceI < nFaces; ++faceI)
                {
                    static string errorMsg =
                        string(" while reading face ")
                            + Foam::name(faceI) + " on patch "
                            + Foam::name(patchI)
                            + " from file " + fName;

                    cueToOrDie(is, "SURF", errorMsg);
                    cueToOrDie(is, "mat", errorMsg);
                    args = cueToOrDie(is, "refs", errorMsg);

                    label nVert = parse<int>(args);

                    face verts(nVert);
                    forAll(verts, vertI)
                    {
                        is.getLine(line);
                        verts[vertI] = parse<int>(line) + patchVertOffset;
                    }

                    if (triangulate && verts.size() > 3)
                    {
                        triFace fTri;

                        // simple face triangulation about f[0].
                        // cannot use face::triangulation
                        // since points are incomplete
                        fTri[0] = verts[0];
                        for (label fp1 = 1; fp1 < verts.size() - 1; ++fp1)
                        {
                            label fp2 = (fp1 + 1) % verts.size();

                            fTri[1] = verts[fp1];
                            fTri[2] = verts[fp2];

                            faceLst.append(fTri);
                            regionLst.append(patchI);
                        }
                    }
                    else
                    {
                        faceLst.append(verts);
                        regionLst.append(patchI);
                    }
                }

                // Done the current patch.
                // Increment the offset vertices are stored at
                patchVertOffset += nPatchPoints;
            }
            else if (cmd == "kids")
            {
                // 'kids' denotes the end of the current patch.
                label nKids = parse<int>(args);

                if (nKids != 0)
                {
                    FatalErrorIn
                    (
                        "fileFormats::AC3DfileFormat::AC3DfileFormat"
                        "(const fileName&)"
                    )
                        << "Can only read objects without kids."
                        << " Encountered " << nKids << " kids when"
                        << " reading patch " << patchI
                        << exit(FatalError);
                }

                // Done reading current patch
                regionNames.insert(patchI, patchName);
                break;
            }
        }
    }

    // transfer to normal lists
    points().transfer(pointLst);
    faces().transfer(faceLst);
    regions().transfer(regionLst);

    setPatches(regionNames);
    stitchFaces(SMALL);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fileFormats::AC3DfileFormat::writeHeader
(
    Ostream& os,
    const List<surfGroup>& patchLst
)
{
    // Write with patches as separate objects under "world" object.
    // Header is taken over from sample file.
    // Defines separate materials for all patches. Recycle colours.

    // Define 8 standard colours as r,g,b components
    static scalar colourMap[] =
    {
        1, 1, 1,
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
        1, 1, 0,
        0, 1, 1,
        1, 0, 1,
        0.5, 0.5, 1
    };

    // Write header. Define materials.
    os  << "AC3Db" << nl;

    forAll(patchLst, patchI)
    {
        const word& pName = patchLst[patchI].name();

        label colourI = patchI % 8;
        label colourCompI = 3 * colourI;

        os  << "MATERIAL \"" << pName << "Mat\" rgb "
            << colourMap[colourCompI] << ' ' << colourMap[colourCompI+1]
            << ' ' << colourMap[colourCompI+2]
            << "  amb 0.2 0.2 0.2  emis 0 0 0  spec 0.5 0.5 0.5  shi 10"
            << "  trans 0"
            << nl;
    }

    os  << "OBJECT world" << nl
        << "kids " << patchLst.size() << endl;
}


void Foam::fileFormats::AC3DfileFormat::write
(
    Ostream& os,
    const keyedSurface& surf
)
{
    labelList faceMap;
    List<surfGroup> patchLst = surf.sortedRegions(faceMap);

    writeHeader(os, patchLst);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        const surfGroup& p = patchLst[patchI];

        os  << "OBJECT poly" << nl
            << "name \"" << p.name() << '"' << endl;

        // Create patch with only patch faces included for ease of addressing
        boolList include(surf.size(), false);

        forAll(p, patchFaceI)
        {
            const label faceI = faceMap[faceIndex++];

            include[faceI] = true;
        }

        labelList pMap;
        labelList fMap;

        keyedSurface patch = surf.subsetMesh
        (
            include, pMap, fMap
        );

        // Now we have triSurface for this patch alone. Write it.
        os << "numvert " << patch.nPoints() << endl;

        forAll(patch.localPoints(), ptI)
        {
            const point& pt = patch.localPoints()[ptI];

            os << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
        }

        os << "numsurf " << patch.localFaces().size() << endl;

        forAll(patch.localFaces(), faceI)
        {
            const face& f = patch.localFaces()[faceI];

            os  << "SURF 0x20" << nl          // polygon
                << "mat " << patchI << nl
                << "refs " << f.size() << nl;

            forAll(f, fp)
            {
                os << f[fp] << " 0 0" << nl;
            }
        }

        os << "kids 0" << endl;
    }
}


void Foam::fileFormats::AC3DfileFormat::write
(
    Ostream& os,
    const meshedSurface& surf
)
{
    const pointField& pointLst = surf.points();
    const List<FaceType>& faceLst = surf.faces();
    const List<surfGroup>& patchLst = surf.patches();

    writeHeader(os, patchLst);

    forAll(patchLst, patchI)
    {
        const surfGroup& p = patchLst[patchI];

        os  << "OBJECT poly" << nl
            << "name \"" << p.name() << '"' << endl;

        // Temporary primitivePatch to calculate compact points & faces
        primitivePatch patch
        (
            SubList<face>
            (
                faceLst,
                p.start(),
                p.size()
            ),
            pointLst
        );

        os << "numvert " << patch.nPoints() << endl;

        forAll(patch.localPoints(), ptI)
        {
            const point& pt = patch.localPoints()[ptI];

            os << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
        }

        os << "numsurf " << patch.localFaces().size() << endl;

        forAll(patch.localFaces(), faceI)
        {
            const face& f = patch.localFaces()[faceI];

            os  << "SURF 0x20" << nl          // polygon
                << "mat " << patchI << nl
                << "refs " << f.size() << nl;

            forAll(f, fp)
            {
                os << f[fp] << " 0 0" << nl;
            }
        }

        os << "kids 0" << endl;
    }
}

// ************************************************************************* //
