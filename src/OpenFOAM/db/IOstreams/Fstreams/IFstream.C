/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IFstream, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IFstream::IFstream
(
    const fileName& pathname,
    IOstreamOption streamOpt
)
:
    Foam::ifstreamPointer(pathname),
    ISstream(*(ifstreamPointer::get()), pathname, streamOpt)
{
    IOstreamOption::compression(ifstreamPointer::whichCompression());

    setClosed();

    setState(ifstreamPointer::get()->rdstate());

    if (good())
    {
        setOpened();
    }
    else
    {
        setBad();
    }

    lineNumber_ = 1;

    if (debug)
    {
        if (pathname.empty())
        {
            InfoInFunction
                << "Cannot open empty file name"
                << Foam::endl;
        }
        else if (IOstreamOption::COMPRESSED == IOstreamOption::compression())
        {
            InfoInFunction
                << "Decompressing " << (this->name() + ".gz") << Foam::endl;
        }

        if (!opened())
        {
            InfoInFunction
                << "Could not open file " << pathname
                << " for input\n" << info() << Foam::endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::istream& Foam::IFstream::stdStream()
{
    std::istream* ptr = ifstreamPointer::get();

    if (!ptr)
    {
        FatalErrorInFunction
            << "No stream allocated\n"
            << abort(FatalError);
    }

    return *ptr;
}


const std::istream& Foam::IFstream::stdStream() const
{
    const std::istream* ptr = ifstreamPointer::get();

    if (!ptr)
    {
        FatalErrorInFunction
            << "No stream allocated\n"
            << abort(FatalError);
    }

    return *ptr;
}


void Foam::IFstream::rewind()
{
    if (IOstreamOption::COMPRESSED == ifstreamPointer::whichCompression())
    {
        lineNumber_ = 1;  // Reset line number
        ifstreamPointer::reopen_gz(this->name() + ".gz");
        setState(ifstreamPointer::get()->rdstate());
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
