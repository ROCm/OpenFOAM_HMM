/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "lumpedPointIOMovement.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lumpedPointIOMovement, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::lumpedPointIOMovement*
Foam::lumpedPointIOMovement::lookupInRegistry(const objectRegistry& obr)
{
    return obr.lookupObjectPtr<lumpedPointIOMovement>
    (
        lumpedPointMovement::dictionaryName
    );
}


Foam::autoPtr<Foam::lumpedPointIOMovement>
Foam::lumpedPointIOMovement::New
(
    const objectRegistry& obr,
    label ownerId
)
{
    return autoPtr<lumpedPointIOMovement>
    (
        new lumpedPointIOMovement
        (
            IOobject
            (
                lumpedPointMovement::dictionaryName,
                obr.time().caseSystem(),
                obr,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                true // register object
            ),
            ownerId  // tag this patch as owner too
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedPointIOMovement::lumpedPointIOMovement
(
    const IOobject& io,
    label ownerId
)
:
    lumpedPointMovement(),
    regIOobject(io)
{
    bool ok =
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
    );

    if (ok)
    {
        ok = readData(readStream(typeName));
        close();

        if (ok)
        {
            this->ownerId(ownerId);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lumpedPointIOMovement::~lumpedPointIOMovement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::lumpedPointIOMovement::readData(Istream& is)
{
    dictionary dict(is);

    readDict(dict);

    return is.check(FUNCTION_NAME);
}


bool Foam::lumpedPointIOMovement::writeData(Ostream& os) const
{
    os  << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const lumpedPointIOMovement& obj)
{
    obj.writeDict(os);

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
