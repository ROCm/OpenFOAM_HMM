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

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "IFstream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Read surf grouping, points, faces directly from Istream
template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::read(Istream& is)
{
    clear();

    // triangulation required?
    bool mustTriangulate = false;
    {
        Face f;
        if (f.max_size() == 3)
        {
            mustTriangulate = true;
        }
    }

    List<surfGroup> patchLst(is);
    // read points:
    is >> storedPoints();


    // read faces:
    // TODO - specialization to triangulate on-the-fly
    if (mustTriangulate)
    {
        is >> storedFaces();
    }
    else
    {
        is >> storedFaces();
    }

    patches_.setSize(patchLst.size());
    regions_.setSize(size());

    // copy patch info and set regions:
    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        patches_[patchI] = PatchRegionType
        (
            patchLst[patchI],
            patchI
        );

        forAll(patchLst[patchI], patchFaceI)
        {
            regions_[faceIndex++] = patchI;
        }
    }

    return is.good();
}


// write (sorted)
template<class Face>
void Foam::UnsortedMeshedSurface<Face>::write(Ostream& os) const
{
    const List<Face>& faceLst = faces();
    labelList faceMap;
    List<surfGroup> patchLst = sortedRegions(faceMap);

    // just emit some information until we get a nice IOobject
    IOobject::writeBanner(os);
    os  << "// OpenFOAM Surface format" << nl
        << "// ~~~~~~~~~~~~~~~~~~~~~~~" << nl
        << "// regions:" << nl
        << patchLst.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    forAll(patchLst, patchI)
    {
        patchLst[patchI].writeDict(os);
    }
    os  << decrIndent << token::END_LIST << nl;

    IOobject::writeDivider(os);

    // Note: Write with global point numbering
    os  << "\n// points:" << nl << points() << nl;

    IOobject::writeDivider(os);
    os  << "\n// faces:"  << nl;
    os  << faceLst.size() << nl << token::BEGIN_LIST << nl;

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        // Print all faces belonging to this region
        const surfGroup& patch = patchLst[patchI];

        forAll(patch, patchFaceI)
        {
            const face& f = faceLst[faceMap[faceIndex++]];

            os  << f << nl;
        }
    }
    os << token::END_LIST << nl;

    IOobject::writeDivider(os);

    // Check state of Ostream
    os.check("meshedSurface::write(Ostream&)");
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Face>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const UnsortedMeshedSurface<Face>& surf
)
{
    surf.write(os);
    return os;
}

// ************************************************************************* //
