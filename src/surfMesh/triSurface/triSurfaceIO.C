/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "triSurface.H"
#include "Fstream.H"
#include "Time.H"
#include "boundBox.H"
#include "PackedBoolList.H"
#include "surfZoneList.H"
#include "surfaceFormatsCore.H"
#include "MeshedSurfaceProxy.H"
#include "MeshedSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::wordHashSet Foam::triSurface::readTypes_;
Foam::wordHashSet Foam::triSurface::writeTypes_;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::wordHashSet& Foam::triSurface::readTypes()
{
    // Stop-gap measure until reading is handled more generally
    if (readTypes_.empty())
    {
        readTypes_ = { "ftr", "stl", "stlb" };
        readTypes_ += UnsortedMeshedSurface<labelledTri>::readTypes();
        readTypes_ += MeshedSurface<labelledTri>::readTypes();
    }

    return readTypes_;
}


const Foam::wordHashSet& Foam::triSurface::writeTypes()
{
    // Stop-gap measure until writing is handled more generally
    if (writeTypes_.empty())
    {
        writeTypes_ = { "ftr", "stl", "stlb", "gts" };
        writeTypes_ += MeshedSurfaceProxy<labelledTri>::writeTypes();
    }

    return writeTypes_;
}


bool Foam::triSurface::canReadType(const word& ext, const bool verbose)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        readTypes(),
        ext,
        verbose,
        "reading"
    );
}


bool Foam::triSurface::canWriteType(const word& ext, const bool verbose)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        writeTypes(),
        ext,
        verbose,
        "writing"
    );
}


bool Foam::triSurface::canRead(const fileName& name, const bool verbose)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }
    return canReadType(ext, verbose);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::triSurface::printTriangle
(
    Ostream& os,
    const string& pre,
    const labelledTri& f,
    const pointField& points
)
{
    os
        << pre.c_str() << "vertex numbers:"
        << f[0] << ' ' << f[1] << ' ' << f[2] << endl
        << pre.c_str() << "vertex coords :"
        << points[f[0]] << ' ' << points[f[1]] << ' ' << points[f[2]]
        << pre.c_str() << "region        :" << f.region() << endl
        << endl;
}


bool Foam::triSurface::read(Istream& is)
{
    // Read triangles, points from Istream
    is  >> patches_ >> storedPoints() >> storedFaces();

    return true;
}


bool Foam::triSurface::read
(
    const fileName& name,
    const word& ext,
    const bool check
)
{
    if (check && !exists(name))
    {
        FatalErrorInFunction
            << "Cannnot read " << name << exit(FatalError);
    }

    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();

        // Do not check for existence. Let IFstream do the unzipping.
        return read(unzipName, unzipName.ext(), false);
    }

    // Hard-coded readers
    if (ext == "ftr")
    {
        return read(IFstream(name)());
    }
    else if (ext == "stl")
    {
        return readSTL(name);  // ASCII
    }
    else if (ext == "stlb")
    {
        return readSTL(name, true); // Force BINARY
    }

    // UnsortedMeshedSurface
    {
        using proxyType = UnsortedMeshedSurface<labelledTri>;
        if (proxyType::readTypes().found(ext))
        {
            reset(proxyType::New(name, ext)());
            return true;
        }
    }

    // MeshedSurface
    {
        using proxyType = MeshedSurface<labelledTri>;
        if (proxyType::readTypes().found(ext))
        {
            reset(proxyType::New(name, ext)());
            return true;
        }
    }


    FatalErrorInFunction
        << "unknown file extension " << ext
        << " for reading file " << name
        << ". Supported extensions:" << nl
        << "    " << flatOutput(readTypes_.sortedToc()) << nl
        << exit(FatalError);

    return false;
}


void Foam::triSurface::write
(
    const fileName& name,
    const word& ext,
    const bool sort
) const
{
    // Hard-coded readers

    if (ext == "ftr")
    {
        OFstream os(name);
        write(os);
    }
    else if (ext == "stl")
    {
        writeSTLASCII(name, sort);
    }
    else if (ext == "stlb")
    {
        writeSTLBINARY(name);
    }
    else if (ext == "gts")
    {
        writeGTS(name, sort);
    }
    else if (MeshedSurfaceProxy<labelledTri>::canWriteType(ext))
    {
        labelList faceMap;
        List<surfZone> zoneLst = this->sortedZones(faceMap);

        MeshedSurfaceProxy<labelledTri> proxy
        (
            this->points(),
            this->surfFaces(),
            zoneLst,
            faceMap
        );

        proxy.write(name, ext);
    }
    else
    {
        FatalErrorInFunction
            << "unknown file extension " << ext
            << " for writing file " << name
            << ". Supported extensions:" << nl
            << "    " << flatOutput(writeTypes_.sortedToc()) << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurface::triSurface(Istream& is)
:
    ParentType(List<Face>(), pointField()),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    read(is);

    setDefaultPatches();
}


Foam::triSurface::triSurface(const Time& d)
:
    ParentType(List<Face>(), pointField()),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    fileName foamFile(d.caseName() + ".ftr");

    fileName foamPath(d.path()/triSurfInstance(d)/typeName/foamFile);

    IFstream foamStream(foamPath);

    read(foamStream);

    setDefaultPatches();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurface::write
(
    const fileName& name,
    const bool sortByRegion
) const
{
    write(name, name.ext(), sortByRegion);
}


void Foam::triSurface::write(Ostream& os) const
{
    os  << patches() << nl;

    //Note: Write with global point numbering
    os  << points() << nl
        << static_cast<const List<labelledTri>&>(*this) << nl;

    // Check state of Ostream
    os.check(FUNCTION_NAME);
}


void Foam::triSurface::write(const Time& d) const
{
    fileName foamFile(d.caseName() + ".ftr");

    fileName foamPath(d.path()/triSurfInstance(d)/typeName/foamFile);

    OFstream foamStream(foamPath);

    write(foamStream);
}


void Foam::triSurface::writeStats(Ostream& os) const
{
    // Unfortunately nPoints constructs meshPoints() so do compact version
    // ourselves.
    PackedBoolList pointIsUsed(points().size());

    label nPoints = 0;
    boundBox bb(boundBox::invertedBox);
    labelHashSet regionsUsed;

    for (const triSurface::FaceType& f : *this)
    {
        regionsUsed.insert(f.region());

        forAll(f, fp)
        {
            const label pointi = f[fp];
            if (pointIsUsed.set(pointi, 1))
            {
                bb.add(points()[pointi]);
                ++nPoints;
            }
        }
    }

    os  << "Triangles    : " << size()
        << " in " << regionsUsed.size() <<  " region(s)" << nl
        << "Vertices     : " << nPoints << nl
        << "Bounding Box : " << bb << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, triSurface& sm)
{
    sm.clearOut();
    sm.read(is);
    sm.setDefaultPatches();
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const triSurface& sm)
{
    sm.write(os);
    return os;
}


// ************************************************************************* //
