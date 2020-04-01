/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "UnsortedMeshedSurface.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
Foam::autoPtr<Foam::UnsortedMeshedSurface<Face>>
Foam::UnsortedMeshedSurface<Face>::New
(
    const fileName& name,
    const word& fileType,
    bool mandatory
)
{
    DebugInFunction
        << "Construct UnsortedMeshedSurface (" << fileType << ")\n";

    auto cstrIter = fileExtensionConstructorTablePtr_->cfind(fileType);

    if (cstrIter.found())
    {
        return autoPtr<UnsortedMeshedSurface<Face>>(cstrIter()(name));
    }


    // Delegate to friend if possible
    const wordHashSet delegate(MeshReference::readTypes());

    if (delegate.found(fileType))
    {
        // OK, can create indirectly
        auto surf = autoPtr<UnsortedMeshedSurface<Face>>::New();
        surf->transfer(*(MeshReference::New(name, fileType)));

        return surf;
    }
    else if (mandatory)
    {
        FatalErrorInFunction
            << "Unknown surface format " << fileType << nl << nl
            << "Valid types:" << nl
            << flatOutput((delegate | readTypes()).sortedToc()) << nl
            << exit(FatalError);
    }

    // Failed, but was optional
    return nullptr;
}


template<class Face>
Foam::autoPtr<Foam::UnsortedMeshedSurface<Face>>
Foam::UnsortedMeshedSurface<Face>::New(const fileName& name)
{
    word ext(name.ext());
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }

    return New(name, ext);
}


// ************************************************************************* //
