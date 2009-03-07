/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "MeshedSurfaceProxy.H"
#include "Time.H"
#include "surfMesh.H"
#include "MeshedSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Face>
Foam::wordHashSet Foam::MeshedSurfaceProxy<Face>::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
bool Foam::MeshedSurfaceProxy<Face>::canWriteType
(
    const word& ext,
    const bool verbose
)
{
    return checkSupport(writeTypes(), ext, verbose, "writing");
}


template<class Face>
void Foam::MeshedSurfaceProxy<Face>::write
(
    const fileName& name,
    const MeshedSurfaceProxy& surf
)
{
    if (debug)
    {
        Info<< "MeshedSurfaceProxy::write"
            "(const fileName&, const MeshedSurfaceProxy&) : "
            "writing to " << name
            << endl;
    }

    word ext = name.ext();

    typename writefileExtensionMemberFunctionTable::iterator mfIter =
        writefileExtensionMemberFunctionTablePtr_->find(ext);

    if (mfIter == writefileExtensionMemberFunctionTablePtr_->end())
    {
        FatalErrorIn
        (
            "MeshedSurfaceProxy::write(const fileName&)"
        )   << "Unknown file extension " << ext << nl << nl
            << "Valid types are :" << endl
            << writeTypes()
            << exit(FatalError);
    }

    mfIter()(name, surf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurfaceProxy<Face>::MeshedSurfaceProxy
(
    const pointField& pointLst,
    const List<Face>& faceLst,
    const List<surfZone>& zoneLst,
    const List<label>& faceMap
)
:
    points_(pointLst),
    faces_(faceLst),
    zones_(zoneLst),
    faceMap_(faceMap)
{}


template<class Face>
Foam::MeshedSurfaceProxy<Face>::MeshedSurfaceProxy
(
    const MeshedSurface<Face>& surf
)
:
    points_(surf.points()),
    faces_(surf.faces()),
    zones_(surf.surfZones()),
    faceMap_(List<label>())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurfaceProxy<Face>::~MeshedSurfaceProxy()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class Face>
// void Foam::MeshedSurfaceProxy<Face>::write
// (
//     const Time& d,
//     const word& surfName
// ) const
// {
//     write(findMeshFile(d, surfName)());
// }

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
