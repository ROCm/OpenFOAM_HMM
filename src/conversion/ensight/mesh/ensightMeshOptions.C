/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

Foam::ensightMesh::options::options(IOstream::streamFormat format)
:
    format_(format),
    lazy_(false),
    noPatches_(false),
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
    return noPatches_ ? true : !patchPatterns_.valid();
}


bool Foam::ensightMesh::options::usePatches() const
{
    return !noPatches_;
}


bool Foam::ensightMesh::options::useFaceZones() const
{
    return faceZonePatterns_.valid();
}


bool Foam::ensightMesh::options::usePatchSelection() const
{
    return noPatches_ ? false : patchPatterns_.valid();
}


void Foam::ensightMesh::options::reset()
{
    noPatches_ = false;
    patchPatterns_.clear();
    faceZonePatterns_.clear();
}


void Foam::ensightMesh::options::lazy(const bool b)
{
    lazy_ = b;
}


void Foam::ensightMesh::options::noPatches(const bool b)
{
    noPatches_ = b;

    if (noPatches_ && patchPatterns_.valid())
    {
        WarningInFunction
            << " existing patch selection disabled"
            << endl;

        patchPatterns_.clear();
    }
}


void Foam::ensightMesh::options::patchSelection
(
    const wordReList& patterns
)
{
    if (noPatches_)
    {
        WarningInFunction
            << " patch selection specified, but noPatches was already active"
            << endl;
    }
    else
    {
        patchPatterns_.reset(new wordReList(patterns));
    }
}


void Foam::ensightMesh::options::faceZoneSelection
(
    const wordReList& patterns
)
{
    faceZonePatterns_.reset(new wordReList(patterns));
}


const Foam::wordReList& Foam::ensightMesh::options::patchSelection() const
{
    if (usePatchSelection())
    {
        return patchPatterns_();
    }
    else
    {
        return wordReList::null();
    }
}


const Foam::wordReList& Foam::ensightMesh::options::faceZoneSelection() const
{
    if (faceZonePatterns_.valid())
    {
        return faceZonePatterns_();
    }
    else
    {
        return wordReList::null();
    }
}


// ************************************************************************* //
