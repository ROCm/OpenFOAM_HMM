/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "MeshedSurfaceProxy.H"

#include "Time.H"
#include "ListOps.H"
#include "surfMesh.H"
#include "OFstream.H"
#include "faceTraits.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
Foam::wordHashSet Foam::MeshedSurfaceProxy<Face>::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


template<class Face>
bool Foam::MeshedSurfaceProxy<Face>::canWriteType
(
    const word& ext,
    const bool verbose
)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        writeTypes(), ext, verbose, "writing"
    );
}


template<class Face>
void Foam::MeshedSurfaceProxy<Face>::write
(
    const fileName& name,
    const MeshedSurfaceProxy& surf,
    const dictionary& options
)
{
    write(name, name.ext(), surf, options);
}


template<class Face>
void Foam::MeshedSurfaceProxy<Face>::write
(
    const fileName& name,
    const word& ext,
    const MeshedSurfaceProxy& surf,
    const dictionary& options
)
{
    if (debug)
    {
        InfoInFunction << "Writing to " << name << endl;
    }

    auto mfIter = writefileExtensionMemberFunctionTablePtr_->cfind(ext);

    if (!mfIter.found())
    {
        FatalErrorInFunction
            << "Unknown file extension " << ext << nl << nl
            << "Valid types:" << nl
            << flatOutput(writeTypes().sortedToc()) << nl
            << exit(FatalError);
    }

    mfIter()(name, surf, options);
}


template<class Face>
void Foam::MeshedSurfaceProxy<Face>::write
(
    const Time& t,
    const word& surfName
) const
{
    // the surface name to be used
    const word name(surfName.size() ? surfName : surfaceRegistry::defaultName);

    if (debug)
    {
        InfoInFunction << "Writing to " << name << endl;
    }


    // The local location
    const fileName objectDir
    (
        t.timePath()/surfaceRegistry::prefix/name/surfMesh::meshSubDir
    );

    if (!isDir(objectDir))
    {
        mkDir(objectDir);
    }


    // write surfMesh/points
    {
        pointIOField io
        (
            IOobject
            (
                "points",
                t.timeName(),
                surfMesh::meshSubDir,
                t,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        OFstream os
        (
            objectDir/io.name(),
            t.writeFormat(),
            IOstream::currentVersion,
            t.writeCompression()
        );

        io.writeHeader(os);

        os  << this->points();

        io.writeEndDivider(os);
    }


    // write surfMesh/faces
    {
        faceCompactIOList io
        (
            IOobject
            (
                "faces",
                t.timeName(),
                surfMesh::meshSubDir,
                t,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        OFstream os
        (
            objectDir/io.name(),
            t.writeFormat(),
            IOstream::currentVersion,
            t.writeCompression()
        );
        io.writeHeader(os);

        if (this->useFaceMap())
        {
            os  << UIndirectList<Face>(this->surfFaces(), this->faceMap());
        }
        else
        {
            os  << this->surfFaces();
        }

        io.writeEndDivider(os);
    }


    // write surfMesh/surfZones
    {
        surfZoneIOList io
        (
            IOobject
            (
                "surfZones",
                t.timeName(),
                surfMesh::meshSubDir,
                t,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        // write as ascii
        OFstream os(objectDir/io.name());
        io.writeHeader(os);

        os  << this->surfZones();

        io.writeEndDivider(os);
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurfaceProxy<Face>::MeshedSurfaceProxy
(
    const pointField& pointLst,
    const UList<Face>& faceLst,
    const UList<surfZone>& zoneLst,
    const labelUList& faceMap
)
:
    points_(pointLst),
    faces_(faceLst),
    zones_(zoneLst),
    faceMap_(faceMap)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
inline Foam::label Foam::MeshedSurfaceProxy<Face>::nTriangles() const
{
    if (faceTraits<Face>::isTri())
    {
        return this->size();
    }

    label nTri = 0;
    for (const Face& f : this->surfFaces())
    {
        nTri += f.nTriangles();
    }

    return nTri;
}


// ************************************************************************* //
