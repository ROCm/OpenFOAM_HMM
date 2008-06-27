/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "solidParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidParticle::solidParticle
(
    const Cloud<solidParticle>& cloud,
    Istream& is,
    bool readFields
)
:
    Particle<solidParticle>(cloud, is)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            d_ = readScalar(is);
            is >> U_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&d_),
                sizeof(d_) + sizeof(U_)
            );
        }
    }

    // Check state of Istream
    is.check("solidParticle::solidParticle(Istream&)");
}


namespace Foam
{

template<>
void Cloud<solidParticle>::readFields()
{
    IOField<scalar> d(fieldIOobject("d"));
    checkFieldIOobject(*this, d);

    IOField<vector> U(fieldIOobject("U"));
    checkFieldIOobject(*this, U);

    label i = 0;
    forAllIter(Cloud<solidParticle>::iterator, *this, iter)
    {
        solidParticle& p = iter();

        p.d_ = d[i];
        p.U_ = U[i];
        i++;
    }
}


template<>
void Cloud<solidParticle>::writeFields() const
{
    label np = size();

    IOField<scalar> d(fieldIOobject("d"), np);
    IOField<vector> U(fieldIOobject("U"), np);

    label i = 0;
    forAllConstIter(Cloud<solidParticle>, *this, iter)
    {
        const solidParticle& p = iter();

        d[i] = p.d_;
        U[i] = p.U_;
        i++;
    }

    d.write();
    U.write();
}

};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const solidParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const Particle<solidParticle>&>(p)
            << token::SPACE << p.d_
            << token::SPACE << p.U_;
    }
    else
    {
        os  << static_cast<const Particle<solidParticle>&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.d_),
            sizeof(p.d_) + sizeof(p.U_)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const solidParticle&)");

    return os;
}


// ************************************************************************* //
