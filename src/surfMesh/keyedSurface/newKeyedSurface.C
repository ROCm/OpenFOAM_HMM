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

#include "keyedSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::autoPtr<Foam::keyedSurface>
Foam::keyedSurface::New
(
    const fileName& fName,
    const word& ext,
    const bool triangulate
)
{
    if (debug)
    {
        Info<< "keyedSurface::New"
            "(const fileName&, const word&, const bool) : "
            "constructing keyedSurface"
            << endl;
    }

    fileExtensionConstructorTable::iterator cstrIter =
        fileExtensionConstructorTablePtr_->find(ext);

    if (cstrIter == fileExtensionConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "keyedSurface::New(const fileName&, const word&, const bool) : "
            "constructing keyedSurface"
        )   << "Unknown file extension " << ext << nl << nl
            << "Valid types are :" << nl
            << fileExtensionConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<keyedSurface>(cstrIter()(fName, triangulate));
}


Foam::autoPtr<Foam::keyedSurface>
Foam::keyedSurface::New
(
    const fileName& fName,
    const bool triangulate
)
{
    const word ext = fName.ext();

    if (debug)
    {
        Info<< "keyedSurface::New(const fileName&, const bool) : "
               "constructing keyedSurface"
            << endl;
    }

    // TODO:
    // if (ext == "gz")
    // {
    //    fileName unzipName = fName.lessExt();
    //
    //    return New(unzipName, unzipName.ext(), ext);
    // }

    return New(fName, ext, triangulate);
}

// ************************************************************************* //
