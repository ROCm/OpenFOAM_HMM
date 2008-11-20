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

#include "AC3DsurfaceFormat.H"
#include "clock.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "tensor.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::AC3DsurfaceFormat<Face>::AC3DsurfaceFormat
(
    const fileName& fName
)
{
    read(fName);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::AC3DsurfaceFormat<Face>::read
(
    const fileName& fName
)
{
    const bool mustTriangulate = this->isTri();
    this->clear();

    IFstream is(fName);
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::AC3DsurfaceFormat::read(const fileName&)"
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
            "fileFormats::AC3DsurfaceFormat::read(const fileName&)"
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
            "fileFormats::AC3DsurfaceFormat::read(const fileName&)"
        )
            << "Cannot find \"OBJECT world\" in file " << fName
            << exit(FatalError);
    }

    // # of kids is the # of patches
    args = cueToOrDie(is, "kids");
    label nPatches = parse<int>(args);

    // Start of vertices for object/patch
    label patchVertOffset = 0;

    DynamicList<point> dynPoints;
    DynamicList<Face>  dynFaces;
    List<label> sizes(nPatches, 0);

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
                    "fileFormats::AC3DsurfaceFormat::read(const fileName&)"
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
                    "fileFormats::AC3DsurfaceFormat::read"
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
                    dynPoints.append(location + pt);
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

                    List<label> verts(nVert);
                    forAll(verts, vertI)
                    {
                        is.getLine(line);
                        verts[vertI] = parse<int>(line) + patchVertOffset;
                    }

                    UList<label>& f = static_cast<UList<label>&>(verts);

                    if (mustTriangulate && f.size() > 3)
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

                            dynFaces.append(fTri);
                            sizes[patchI]++;
                        }
                    }
                    else
                    {
                        dynFaces.append(Face(f));
                        sizes[patchI]++;
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
                        "fileFormats::AC3DsurfaceFormat::read(const fileName&)"
                    )
                        << "Can only read objects without kids."
                        << " Encountered " << nKids << " kids when"
                        << " reading patch " << patchI
                        << exit(FatalError);
                }

                // Done reading current patch
                break;
            }
        }
    }

    // transfer to normal lists
    this->storedPoints().transfer(dynPoints);
    this->storedFaces().transfer(dynFaces);

    this->addPatches(sizes);
    this->stitchFaces(SMALL);
    return true;
}


template<class Face>
void Foam::fileFormats::AC3DsurfaceFormat<Face>::write
(
    Ostream& os,
    const MeshedSurface<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>& faceLst = surf.faces();
    const List<surfGroup>& patchLst = surf.patches();

    writeHeader(os, patchLst);

    forAll(patchLst, patchI)
    {
        const surfGroup& p = patchLst[patchI];

        os  << "OBJECT poly" << nl
            << "name \"" << p.name() << '"' << endl;

        // Temporary PrimitivePatch to calculate compact points & faces
        // use 'UList' to avoid allocations!
        PrimitivePatch<Face, UList, const pointField&> patch
        (
            faceLst,
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
            const Face& f = patch.localFaces()[faceI];

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


template<class Face>
void Foam::fileFormats::AC3DsurfaceFormat<Face>::write
(
    Ostream& os,
    const UnsortedMeshedSurface<Face>& surf
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

        UnsortedMeshedSurface<Face> subm = surf.subsetMesh(include);

        // Now we have isolated surface for this patch alone. Write it.
        os << "numvert " << subm.nPoints() << endl;

        forAll(subm.localPoints(), ptI)
        {
            const point& pt = subm.localPoints()[ptI];

            os << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
        }

        os << "numsurf " << subm.localFaces().size() << endl;

        forAll(subm.localFaces(), faceI)
        {
            const Face& f = subm.localFaces()[faceI];

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
