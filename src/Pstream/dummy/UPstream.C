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

#include "UPstream.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{}


bool Foam::UPstream::init(int& argc, char**& argv)
{
    FatalErrorIn("UPstream::init(int& argc, char**& argv)")
        << "Trying to use the dummy Pstream library." << nl
        << "This dummy library cannot be used in parallel mode"
        << Foam::exit(FatalError);

    return false;
}


void Foam::UPstream::exit(int errnum)
{
    notImplemented("UPstream::exit(int errnum)");
}


void Foam::UPstream::abort()
{
    notImplemented("UPstream::abort()");
}


void Foam::reduce(scalar&, const sumOp<scalar>&)
{}



void Foam::UPstream::waitRequests()
{}


bool Foam::UPstream::finishedRequest(const label i)
{
    notImplemented("UPstream::finishedRequest()");
    return false;
}


// ************************************************************************* //
