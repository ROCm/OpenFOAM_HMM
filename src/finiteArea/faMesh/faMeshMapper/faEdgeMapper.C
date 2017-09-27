/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

Description
    FA edge mapper.

\*---------------------------------------------------------------------------*/

#include "faEdgeMapper.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faEdgeMapper::calcAddressing() const
{
    if (directAddrPtr_)
    {
        FatalErrorInFunction
            << "Addressing already calculated"
            << abort(FatalError);
    }

    hasUnmapped_ = false;

    // Dummy mapping: take value from edge 0
    directAddrPtr_ = new labelList(size(), 0);
}


void Foam::faEdgeMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    hasUnmapped_ = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faEdgeMapper::faEdgeMapper
(
    const faMesh& mesh,
    const mapPolyMesh& mpm
)
:
    mesh_(mesh),
    mpm_(mpm),
    sizeBeforeMapping_(mesh.nInternalEdges()),
    hasUnmapped_(false),
    directAddrPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faEdgeMapper::~faEdgeMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::faEdgeMapper::directAddressing() const
{
    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList& Foam::faEdgeMapper::addressing() const
{
    FatalErrorInFunction
        << "Requested interpolative addressing for a direct mapper."
        << abort(FatalError);

    return labelListList::null();
}


const Foam::scalarListList& Foam::faEdgeMapper::weights() const
{
    FatalErrorInFunction
        << "Requested interpolative weights for a direct mapper."
        << abort(FatalError);

    return scalarListList::null();
}


// ************************************************************************* //
