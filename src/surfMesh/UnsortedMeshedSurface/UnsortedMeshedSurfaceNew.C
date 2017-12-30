/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    const word& ext
)
{
    if (debug)
    {
        InfoInFunction << "Constructing UnsortedMeshedSurface" << endl;
    }

    auto cstrIter = fileExtensionConstructorTablePtr_->cfind(ext);

    if (!cstrIter.found())
    {
        // No direct reader, delegate to parent if possible
        const wordHashSet& delegate = ParentType::readTypes();

        if (delegate.found(ext))
        {
            // Create indirectly
            autoPtr<UnsortedMeshedSurface<Face>> surf
            (
                new UnsortedMeshedSurface<Face>
            );
            surf().transfer(ParentType::New(name, ext)());

            return surf;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown file extension " << ext << nl << nl
                << "Valid types:" << nl
                << flatOutput((delegate | readTypes()).sortedToc()) << nl
                << exit(FatalError);
        }
    }

    return autoPtr<UnsortedMeshedSurface<Face>>(cstrIter()(name));
}


template<class Face>
Foam::autoPtr<Foam::UnsortedMeshedSurface<Face>>
Foam::UnsortedMeshedSurface<Face>::New(const fileName& name)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }

    return New(name, ext);
}


// ************************************************************************* //
