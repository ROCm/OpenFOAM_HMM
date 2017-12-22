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
    // Note: does not include position_
    sizeof(label) + sizeof(scalar) + sizeof(scalar) + sizeof(vector)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::injectedParticle::injectedParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, false), // force to read old positions file
    position_(Zero),
    tag_(-1),
    soi_(0.0),
    d_(0.0),
    U_(Zero)
{
    if (readFields)
    {
        // After the base particle class has read the fields from file and
        // constructed the necessary barcentric co-ordinates we can update the
        // particle position on this mesh
        position_ = particle::position();

        if (is.format() == IOstream::ASCII)
        {
            tag_ = readLabel(is);
            soi_ = readScalar(is);
            d_ = readScalar(is);
            is >> U_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&soi_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


void Foam::injectedParticle::readFields(Cloud<injectedParticle>& c)
{
    if (!c.size())
    {
        return;
    }

    // Note: not reading local position_ - defer to base particle class

    particle::readFields(c);

    IOField<label> tag(c.fieldIOobject("tag", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tag);

    IOField<scalar> soi(c.fieldIOobject("soi", IOobject::MUST_READ));
    c.checkFieldIOobject(c, soi);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    label i = 0;

    forAllIters(c, iter)
    {
        injectedParticle& p = iter();

        p.tag_ = tag[i];
        p.soi_ = soi[i];
        p.d_ = d[i];
        p.U_ = U[i];

        i++;
    }
}


void Foam::injectedParticle::writeFields(const Cloud<injectedParticle>& c)
{
    // Force writing positions instead of coordinates
    particle::writeLagrangianCoordinates = false;
    particle::writeLagrangianPositions = true;

    particle::writeFields(c);

    // Note: not writing local position_ - defer to base particle class

    label np =  c.size();

    IOField<label> tag(c.fieldIOobject("tag", IOobject::NO_READ), np);
    IOField<scalar> soi(c.fieldIOobject("soi", IOobject::NO_READ), np);
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIters(c, iter)
    {
        const injectedParticle& p = iter();

        tag[i] = p.tag();
        soi[i] = p.soi();
        d[i] = p.d();
        U[i] = p.U();

        i++;
    }

    tag.write();
    soi.write();
    d.write();
    U.write();
}


void Foam::injectedParticle::writeObjects
(
    const Cloud<injectedParticle>& c,
    objectRegistry& obr
)
{
    // Force writing positions instead of coordinates
    particle::writeLagrangianCoordinates = false;
    particle::writeLagrangianPositions = true;

    particle::writeObjects(c, obr);

    label np = c.size();

    IOField<label>& tag(cloud::createIOField<label>("tag", np, obr));
    IOField<scalar>& soi(cloud::createIOField<scalar>("soi", np, obr));
    IOField<scalar>& d(cloud::createIOField<scalar>("d", np, obr));
    IOField<vector>& U(cloud::createIOField<vector>("U", np, obr));

    label i = 0;

    forAllConstIters(c, iter)
    {
        const injectedParticle& p = iter();

        tag[i] = p.tag();
        soi[i] = p.soi();
        d[i] = p.d();
        U[i] = p.U();

        i++;
    }
}


void Foam::injectedParticle::writePosition(Ostream& os) const
{
    if (os.format() == IOstream::ASCII)
    {
        os  << position_ << token::SPACE << cell();
    }
    else
    {
        struct oldParticle
        {
            vector position;
            label celli;
            label facei;
            scalar stepFraction;
            label tetFacei;
            label tetPti;
            label origProc;
            label origId;
        } p;

        const size_t s =
            offsetof(oldParticle, facei) - offsetof(oldParticle, position);

        p.position = position_;
        p.celli = cell();

        os.write(reinterpret_cast<const char*>(&p.position), s);
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const injectedParticle& p
)
{
    // Note: not writing local position_ - defer to base particle class

    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.tag()
            << token::SPACE << p.soi()
            << token::SPACE << p.d()
            << token::SPACE << p.U();
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.tag_),
            injectedParticle::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
