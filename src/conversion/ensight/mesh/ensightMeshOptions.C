/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "ensightMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightMesh::options::options()
:
    options(IOstream::streamFormat::BINARY)
{}


Foam::ensightMesh::options::options(IOstream::streamFormat format)
:
    format_(format),
    lazy_(false),
    internal_(true),
    boundary_(true),
    patchPatterns_(),
    faceZonePatterns_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::IOstream::streamFormat Foam::ensightMesh::options::format() const
{
    return format_;
}


bool Foam::ensightMesh::options::lazy() const
{
    return lazy_;
}


bool Foam::ensightMesh::options::useInternalMesh() const
{
    return internal_;
}


bool Foam::ensightMesh::options::useBoundaryMesh() const
{
    return boundary_;
}


bool Foam::ensightMesh::options::useFaceZones() const
{
    return faceZonePatterns_.size();
}


void Foam::ensightMesh::options::reset()
{
    internal_ = true;
    boundary_ = true;
    patchPatterns_.clear();
    faceZonePatterns_.clear();
}


void Foam::ensightMesh::options::lazy(bool beLazy)
{
    lazy_ = beLazy;
}


void Foam::ensightMesh::options::useInternalMesh(bool on)
{
    internal_ = on;
}


void Foam::ensightMesh::options::useBoundaryMesh(bool on)
{
    boundary_ = on;

    if (!boundary_ && patchPatterns_.size())
    {
        patchPatterns_.clear();

        WarningInFunction
            << "Deactivating boundary and removing old patch selection"
            << endl;
    }
}


void Foam::ensightMesh::options::patchSelection
(
    const UList<wordRe>& patterns
)
{
    patchPatterns_ = wordRes(patterns);

    if (!boundary_ && patchPatterns_.size())
    {
        patchPatterns_.clear();

        WarningInFunction
            << "Ignoring patch selection, boundary is not active"
            << endl;
    }
}


void Foam::ensightMesh::options::patchSelection
(
    List<wordRe>&& patterns
)
{
    patchPatterns_ = wordRes(std::move(patterns));

    if (!boundary_ && patchPatterns_.size())
    {
        patchPatterns_.clear();

        WarningInFunction
            << "Ignoring patch selection, boundary is not active"
            << endl;
    }
}


void Foam::ensightMesh::options::faceZoneSelection
(
    const UList<wordRe>& patterns
)
{
    faceZonePatterns_ = wordRes(patterns);
}


void Foam::ensightMesh::options::faceZoneSelection
(
    List<wordRe>&& patterns
)
{
    faceZonePatterns_ = wordRes(std::move(patterns));
}


const Foam::wordRes& Foam::ensightMesh::options::patchSelection() const
{
    return patchPatterns_;
}


const Foam::wordRes& Foam::ensightMesh::options::faceZoneSelection() const
{
    return faceZonePatterns_;
}


// ************************************************************************* //
