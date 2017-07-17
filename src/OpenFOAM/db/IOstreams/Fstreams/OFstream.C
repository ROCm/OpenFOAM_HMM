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

#include "OFstream.H"
#include "OSspecific.H"
#include "gzstream.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OFstream, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::OFstreamAllocator::OFstreamAllocator
(
    const fileName& pathname,
    IOstream::compressionType compression
)
:
    allocatedPtr_(nullptr)
{
    if (pathname.empty())
    {
        if (OFstream::debug)
        {
            InfoInFunction << "Cannot open null file " << endl;
        }
    }

    if (compression == IOstream::COMPRESSED)
    {
        // Get identically named uncompressed version out of the way
        if (isFile(pathname, false))
        {
            rm(pathname);
        }

        allocatedPtr_ = new ogzstream((pathname + ".gz").c_str());
    }
    else
    {
        // Get identically named compressed version out of the way
        if (isFile(pathname + ".gz", false))
        {
            rm(pathname + ".gz");
        }

        allocatedPtr_ = new std::ofstream(pathname);
    }
}


Foam::OFstreamAllocator::~OFstreamAllocator()
{
    deallocate();
}


void Foam::OFstreamAllocator::deallocate()
{
    if (allocatedPtr_)
    {
        delete allocatedPtr_;
        allocatedPtr_ = nullptr;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstream::OFstream
(
    const fileName& pathname,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OFstreamAllocator(pathname, compression),
    OSstream(*allocatedPtr_, pathname, format, version, compression)
{
    setClosed();
    setState(allocatedPtr_->rdstate());

    if (!good())
    {
        if (debug)
        {
            InfoInFunction
                << "Could not open file " << pathname
                << " for output" << nl << info() << Foam::endl;
        }

        setBad();
    }
    else
    {
        setOpened();
    }

    lineNumber_ = 1;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OFstream::~OFstream()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::ostream& Foam::OFstream::stdStream()
{
    if (!allocatedPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated." << abort(FatalError);
    }
    return *allocatedPtr_;
}


const std::ostream& Foam::OFstream::stdStream() const
{
    if (!allocatedPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated." << abort(FatalError);
    }
    return *allocatedPtr_;
}


void Foam::OFstream::print(Ostream& os) const
{
    os  << "OFstream: ";
    OSstream::print(os);
}


// ************************************************************************* //
