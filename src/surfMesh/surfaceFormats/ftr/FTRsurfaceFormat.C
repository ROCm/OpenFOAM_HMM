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

#include "FTRsurfaceFormat.H"
#include "Keyed.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::FTRsurfaceFormat<Face>::FTRsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::FTRsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    this->clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::FTRsurfaceFormat::read(const fileName&)"
        )
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    List<ftrPatch> readPatches(is);

    // read points directly
    is >> this->storedPoints();

    // read faces with keys
    List<Keyed<triFace> > readFaces(is);

    List<Face>  faceLst(readFaces.size());
    List<label> regionLst(readFaces.size());

    // disentangle faces/keys - already triangulated
    forAll(readFaces, faceI)
    {
        // unfortunately cannot transfer to save memory
        faceLst[faceI]   = readFaces[faceI];
        regionLst[faceI] = readFaces[faceI].key();
    }

    this->storedFaces().transfer(faceLst);
    this->storedRegions().transfer(regionLst);

    // cast ftrPatch into new form
    List<surfPatchIdentifier> newPatches(readPatches.size());
    forAll(newPatches, patchI)
    {
        newPatches[patchI] = surfPatchIdentifier
        (
            readPatches[patchI].name(),
            patchI
        );
    }

    this->storedPatches().transfer(newPatches);
    return true;
}

// ************************************************************************* //
