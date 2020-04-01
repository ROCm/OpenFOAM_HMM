/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "edgeMesh.H"
#include "boundBox.H"
#include "edgeMeshFormat.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeMesh::edgeMesh
(
    const fileName& name,
    const word& fileType
)
:
    points_(),
    edges_(),
    pointEdgesPtr_(nullptr)
{
    read(name, fileType);
}


Foam::edgeMesh::edgeMesh(const fileName& name)
:
    points_(),
    edges_(),
    pointEdgesPtr_(nullptr)
{
    read(name);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::edgeMesh::read(const fileName& name)
{
    word ext(name.ext());
    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();
        return read(unzipName, unzipName.ext());
    }

    return read(name, ext);
}


bool Foam::edgeMesh::read
(
    const fileName& name,
    const word& fileType
)
{
    // Read via selector mechanism
    transfer(*New(name, fileType));
    return true;
}


void Foam::edgeMesh::write
(
    const fileName& name,
    const word& fileType,
    const edgeMesh& mesh
)
{
    DebugInFunction << "Writing to " << name << endl;

    auto mfIter = writefileExtensionMemberFunctionTablePtr_->cfind(fileType);

    if (!mfIter.found())
    {
        FatalErrorInLookup
        (
            "extension",
            fileType,
            *writefileExtensionMemberFunctionTablePtr_
        ) << exit(FatalError);
    }

    mfIter()(name, mesh);
}


void Foam::edgeMesh::write
(
    const fileName& name,
    const edgeMesh& mesh
)
{
    write(name, name.ext(), mesh);
}


void Foam::edgeMesh::writeStats(Ostream& os) const
{
    os  << indent << "points      : " << points().size() << nl;
    os  << indent << "edges       : " << edges().size() << nl;
    os  << indent << "boundingBox : " << boundBox(this->points()) << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const edgeMesh& em)
{
    fileFormats::edgeMeshFormat::write(os, em.points_, em.edges_);

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, edgeMesh& em)
{
    fileFormats::edgeMeshFormat::read(is, em.points_, em.edges_);

    em.pointEdgesPtr_.reset(nullptr);

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
