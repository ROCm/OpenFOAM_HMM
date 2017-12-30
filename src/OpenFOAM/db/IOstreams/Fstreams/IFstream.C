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

#include "IFstream.H"
#include "OSspecific.H"
#include "gzstream.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IFstream, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IFstreamAllocator::IFstreamAllocator(const fileName& pathname)
:
    allocatedPtr_(nullptr),
    compression_(IOstream::UNCOMPRESSED)
{
    if (pathname.empty())
    {
        if (IFstream::debug)
        {
            InfoInFunction << "Cannot open null file " << endl;
        }
    }

    allocatedPtr_ = new std::ifstream(pathname);

    // If the file is compressed, decompress it before reading.
    if (!allocatedPtr_->good() && isFile(pathname + ".gz", false))
    {
        if (IFstream::debug)
        {
            InfoInFunction << "Decompressing " << pathname + ".gz" << endl;
        }

        delete allocatedPtr_;

        allocatedPtr_ = new igzstream((pathname + ".gz").c_str());

        if (allocatedPtr_->good())
        {
            compression_ = IOstream::COMPRESSED;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IFstreamAllocator::~IFstreamAllocator()
{
    deallocate();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IFstreamAllocator::deallocate()
{
    if (allocatedPtr_)
    {
        delete allocatedPtr_;
        allocatedPtr_ = nullptr;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IFstream::IFstream
(
    const fileName& pathname,
    streamFormat format,
    versionNumber version
)
:
    IFstreamAllocator(pathname),
    ISstream
    (
        *allocatedPtr_,
        pathname,
        format,
        version,
        IFstreamAllocator::compression_
    )
{
    setClosed();

    setState(allocatedPtr_->rdstate());

    if (!good())
    {
        if (debug)
        {
            InfoInFunction
                << "Could not open file " << pathname
                << " for input" << nl << info() << Foam::endl;
        }

        setBad();
    }
    else
    {
        setOpened();
    }

    lineNumber_ = 1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IFstream::~IFstream()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::istream& Foam::IFstream::stdStream()
{
    if (!allocatedPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated"
            << abort(FatalError);
    }
    return *allocatedPtr_;
}


const std::istream& Foam::IFstream::stdStream() const
{
    if (!allocatedPtr_)
    {
        FatalErrorInFunction
            << "No stream allocated"
            << abort(FatalError);
    }
    return *allocatedPtr_;
}


void Foam::IFstream::rewind()
{
    lineNumber_ = 1;      // Reset line number

    igzstream* gzPtr = nullptr;

    try
    {
        gzPtr = dynamic_cast<igzstream*>(allocatedPtr_);
    }
    catch (std::bad_cast)
    {
        gzPtr = nullptr;
    }

    if (gzPtr)
    {
        // Need special treatment for gzstream.
        gzPtr->close();
        gzPtr->clear();
        gzPtr->open((this->name() + ".gz").c_str());

        setState(gzPtr->rdstate());
    }
    else
    {
        ISstream::rewind();
    }
}


void Foam::IFstream::print(Ostream& os) const
{
    os  << "IFstream: ";
    ISstream::print(os);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

Foam::IFstream& Foam::IFstream::operator()() const
{
    if (!good())
    {
        // Also checks .gz file
        if (isFile(this->name(), true))
        {
            check(FUNCTION_NAME);
            FatalIOError.exit();
        }
        else
        {
            FatalIOErrorInFunction(*this)
                << "file " << this->name() << " does not exist"
                << exit(FatalIOError);
        }
    }

    return const_cast<IFstream&>(*this);
}


// ************************************************************************* //
