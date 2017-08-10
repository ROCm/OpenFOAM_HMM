/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd
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

#include "DTRMParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::DTRMParticle::propertyList_ =
    Foam::DTRMParticle::propertyList();

const std::size_t Foam::DTRMParticle::sizeofFields_
(
    sizeof(DTRMParticle) - sizeof(particle)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DTRMParticle::DTRMParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    p0_(position_),
    p1_(point::zero),
    I0_(0),
    I_(0),
    dA_(0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> p0_ >> p1_ >> I0_ >> I_ >> dA_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&p0_), sizeofFields_);
        }
    }
}


Foam::Ostream& Foam::operator<<(Ostream& os, const DTRMParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.p0_
            << token::SPACE << p.p1_
            << token::SPACE << p.I0_
            << token::SPACE << p.I_
            << token::SPACE << p.dA_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.p0_),
            DTRMParticle::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const DTRMParticle&)");

    return os;
}


// ************************************************************************* //
