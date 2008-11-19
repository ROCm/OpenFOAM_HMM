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

#include "OFFsurfaceFormat.H"
#include "IFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::OFFsurfaceFormat<Face>::OFFsurfaceFormat
(
    const fileName& fName
)
{
    read(fName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::OFFsurfaceFormat<Face>::read
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
            "fileFormats::OFFsurfaceFormat<Face>::read(const fileName&)"
        )
            << "Cannot read file " << fName
            << exit(FatalError);
    }

    // Read header
    string hdr = this->getLineNoComment(is);
    if (hdr != "OFF")
    {
        FatalErrorIn
        (
            "fileFormats::OFFsurfaceFormat<Face>::read(const fileName&)"
        )
            << "OFF file " << fName << " does not start with 'OFF'"
            << exit(FatalError);
    }


    // get dimensions
    label nPoints, nElems, nEdges;

    string line = this->getLineNoComment(is);
    {
        IStringStream lineStream(line);
        lineStream >> nPoints >> nElems >> nEdges;
    }

    // Read points
    pointField pointLst(nPoints);
    forAll(pointLst, pointI)
    {
        scalar x, y, z;
        line = this->getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream >> x >> y >> z;
        }
        pointLst[pointI] = point(x, y, z);
    }

    // Read faces - ignore optional region information
    // use a DynamicList for possible on-the-fly triangulation
    DynamicList<Face>  faceLst(nElems);

    for (label faceI = 0; faceI < nElems; ++faceI)
    {
        line = this->getLineNoComment(is);

        {
            IStringStream lineStream(line);

            label nVerts;
            lineStream >> nVerts;

            List<label> verts(nVerts);

            forAll(verts, vertI)
            {
                lineStream >> verts[vertI];
            }

            UList<label>& f = static_cast<UList<label>&>(verts);

            if (mustTriangulate && f.size() > 3)
            {
                triFace fTri;

                // simple face triangulation about f[0].
                fTri[0] = f[0];
                for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
                {
                    label fp2 = (fp1 + 1) % f.size();

                    fTri[1] = f[fp1];
                    fTri[2] = f[fp2];

                    faceLst.append(fTri);
                }
            }
            else
            {
                faceLst.append(Face(f));
            }
        }
    }

    // transfer to normal lists
    reset
    (
        xferMove(pointLst),
        xferMoveTo<List<Face> >(faceLst)
    );

    // no region information
    this->onePatch();
    return true;
}


template<class Face>
void Foam::fileFormats::OFFsurfaceFormat<Face>::write
(
    Ostream& os,
    const UnsortedMeshedSurface<Face>& surf
)
{
    const List<Face>& faceLst = surf.faces();

    labelList faceMap;
    List<surfGroup> patchLst = surf.sortedRegions(faceMap);

    writeHeader(os, surf.points(), faceLst.size(), patchLst);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        os << "# <patch name=\"" << patchLst[patchI].name() << "\">" << endl;

        forAll(patchLst[patchI], patchFaceI)
        {
            const Face& f = faceLst[faceMap[faceIndex++]];

            os << f.size();
            forAll(f, fp)
            {
                os << ' ' << f[fp];
            }

            // add optional region information
            os << ' ' << patchI << endl;
        }
        os << "# </patch>" << endl;
    }
    os << "# </faces>" << endl;
}


template<class Face>
void Foam::fileFormats::OFFsurfaceFormat<Face>::write
(
    Ostream& os,
    const MeshedSurface<Face>& surf
)
{
    const List<Face>& faceLst = surf.faces();
    const List<surfGroup>& patchLst = surf.patches();

    writeHeader(os, surf.points(), faceLst.size(), patchLst);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        os << "# <patch name=\"" << patchLst[patchI].name() << "\">" << endl;

        forAll(patchLst[patchI], patchFaceI)
        {
            const Face& f = faceLst[faceIndex++];

            os << f.size();
            forAll(f, fp)
            {
                os << ' ' << f[fp];
            }

            // add optional region information
            os << ' ' << patchI << endl;
        }
        os << "# </patch>" << endl;
    }
    os << "# </faces>" << endl;
}

// ************************************************************************* //
