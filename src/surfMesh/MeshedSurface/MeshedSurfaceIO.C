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
bool Foam::MeshedSurface<Face>::read(Istream& is)
{
    clear();

    List<surfGroup> patchLst(is);

    // copy and set the indices
    patches_.setSize(patchLst.size());
    forAll(patchLst, patchI)
    {
        patches_[patchI] = surfGroup
        (
            patchLst[patchI],
            patchI
        );
    }

    // read points:
    is >> this->storedPoints();

    // must triangulate?
    if (this->isTri())
    {
        List<face> faceLst(is);

        MeshedSurface<face> surf;
        surf.reset
        (
            xfer<pointField>::null(),
            xferMove(faceLst)
        );
        surf.addPatches(patches_);

        // this will break if the triangulation needed points
        surf.triangulate();
        patches_ = surf.patches();

        // transcribe from face -> triFace (Face)
        const List<face>& origFaces = surf.faces();
        List<Face>  newFaces(origFaces.size());
        forAll(origFaces, faceI)
        {
            newFaces[faceI] = Face
            (
                static_cast<const UList<label>&>(origFaces[faceI])
            );
        }
        surf.clear();

        this->storedFaces().transfer(newFaces);
    }
    else
    {
        // read faces:
        is >> this->storedFaces();
    }

    return is.good();
}


template<class Face>
void Foam::MeshedSurface<Face>::write(Ostream& os) const
{
    // just emit some information until we get a nice IOobject
    IOobject::writeBanner(os);
    os  << "// OpenFOAM Surface Format" << nl
        << "// ~~~~~~~~~~~~~~~~~~~~~~~" << nl
        << "// regions:" << nl
        << patches_.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    forAll(patches_, patchI)
    {
        patches_[patchI].writeDict(os);
    }
    os  << decrIndent << token::END_LIST << nl;

    IOobject::writeDivider(os);

    // Note: Write with global point numbering
    os  << "\n// points:" << nl << this->points() << nl;

    IOobject::writeDivider(os);
    os  << "\n// faces:"  << nl << this->faces() << nl;

    IOobject::writeDivider(os);

    // Check state of Ostream
    os.check("MeshedSurface::write(Ostream&)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Face>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MeshedSurface<Face>& surf
)
{
    surf.write(os);
    return os;
}

// ************************************************************************* //
