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

#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Face>
Foam::autoPtr<Foam::UnsortedMeshedSurface<Face> >
Foam::UnsortedMeshedSurface<Face>::New
(
    const fileName& fName,
    const word& ext
)
{
    if (debug)
    {
        Info<< "UnsortedMeshedSurface<Face>::New"
            "(const fileName&, const word&) : "
            "constructing UnsortedMeshedSurface"
            << endl;
    }

    typename fileExtensionConstructorTable::iterator cstrIter =
        fileExtensionConstructorTablePtr_->find(ext);

    if (cstrIter == fileExtensionConstructorTablePtr_->end())
    {
        // no direct reader, delegate if possible
        wordHashSet supported = SiblingType::readTypes();
        if (supported.found(ext))
        {
            // create indirectly
            autoPtr<UnsortedMeshedSurface<Face> > surf
            (
                new UnsortedMeshedSurface<Face>
            );
            surf().transfer(SiblingType::New(fName, ext)());

            return surf;
        }

        // nothing left but to issue an error
        supported += readTypes();
        supported.insert(nativeExt);

        FatalErrorIn
        (
            "UnsortedMeshedSurface<Face>::New"
            "(const fileName&, const word&) : "
            "constructing UnsortedMeshedSurface"
        )   << "Unknown file extension " << ext << nl << nl
            << "Valid types are :" << nl
            << supported
            << exit(FatalError);
    }

    return autoPtr<UnsortedMeshedSurface<Face> >(cstrIter()(fName));
}


template<class Face>
Foam::autoPtr<Foam::UnsortedMeshedSurface<Face> >
Foam::UnsortedMeshedSurface<Face>::New
(
    const fileName& fName
)
{
    word ext = fName.ext();
    if (ext == "gz")
    {
        ext = fName.lessExt().ext();
    }

    return New(fName, ext);
}

// ************************************************************************* //
