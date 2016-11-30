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

#include "injectedParticle.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::injectedParticle::propertyList_ =
    Foam::injectedParticle::propertyList();

Foam::string Foam::injectedParticle::propertyTypes_ =
    Foam::injectedParticle::propertyTypes();

const std::size_t Foam::injectedParticle::sizeofFields
(
    sizeof(scalar) + sizeof(scalar) + sizeof(vector)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::injectedParticle::injectedParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    soi_(0.0),
    d_(0.0),
    U_(Zero)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            soi_ = readScalar(is);
            d_ = readScalar(is);
            is >> U_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&soi_), sizeofFields);
        }
    }

    // Check state of Istream
    is.check
    (
        "injectedParticle::injectedParticle"
        "(const polyMesh&, Istream&, bool)"
    );
}


void Foam::injectedParticle::readFields(Cloud<injectedParticle>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<scalar> soi(c.fieldIOobject("soi", IOobject::MUST_READ));
    c.checkFieldIOobject(c, soi);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    label i = 0;

    forAllIter(Cloud<injectedParticle>, c, iter)
    {
        injectedParticle& p = iter();

        p.soi_ = soi[i];
        p.d_ = d[i];
        p.U_ = U[i];

        i++;
    }
}


void Foam::injectedParticle::writeFields(const Cloud<injectedParticle>& c)
{
    particle::writeFields(c);

    label np =  c.size();

    IOField<scalar> soi(c.fieldIOobject("soi", IOobject::NO_READ), np);
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(Cloud<injectedParticle>, c, iter)
    {
        const injectedParticle& p = iter();

        soi[i] = p.soi();
        d[i] = p.d();
        U[i] = p.U();

        i++;
    }

    soi.write();
    d.write();
    U.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const injectedParticle& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.soi()
            << token::SPACE << p.d()
            << token::SPACE << p.U();
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.soi_),
            injectedParticle::sizeofFields
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const injectedParticle&)"
    );

    return os;
}


// ************************************************************************* //
