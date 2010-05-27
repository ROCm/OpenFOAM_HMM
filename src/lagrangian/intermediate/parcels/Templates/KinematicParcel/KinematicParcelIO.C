/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "KinematicParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class ParcelType>
Foam::string Foam::KinematicParcel<ParcelType>::propHeader =
    Particle<ParcelType>::propHeader
  + " active"
  + " typeId"
  + " nParticle"
  + " d"
  + " (Ux Uy Uz)"
  + " (fx fy fz)"
  + " (angularMomentumx angularMomentumy angularMomentumz)"
  + " (torquex torquey torquez)"
  + " rho"
  + " tTurb"
  + " (UTurbx UTurby UTurbz)"
  + " pairCollisions wallCollisions";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const Cloud<ParcelType>& cloud,
    Istream& is,
    bool readFields
)
:
    Particle<ParcelType>(cloud, is, readFields),
    active_(false),
    typeId_(0),
    nParticle_(0.0),
    d_(0.0),
    U_(vector::zero),
    f_(vector::zero),
    angularMomentum_(vector::zero),
    torque_(vector::zero),
    rho_(0.0),
    tTurb_(0.0),
    UTurb_(vector::zero),
    collisionRecords_(),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            active_ = readBool(is);
            typeId_ = readLabel(is);
            nParticle_ = readScalar(is);
            d_ = readScalar(is);
            is >> U_;
            is >> f_;
            is >> angularMomentum_;
            is >> torque_;
            rho_ = readScalar(is);
            tTurb_ = readScalar(is);
            is >> UTurb_;
            is >> collisionRecords_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&active_),
                sizeof(active_)
              + sizeof(typeId_)
              + sizeof(nParticle_)
              + sizeof(d_)
              + sizeof(U_)
              + sizeof(f_)
              + sizeof(angularMomentum_)
              + sizeof(torque_)
              + sizeof(rho_)
              + sizeof(tTurb_)
              + sizeof(UTurb_)
            );
            is  >> collisionRecords_;
        }
    }

    // Check state of Istream
    is.check
    (
        "KinematicParcel<ParcelType>::KinematicParcel"
        "(const Cloud<ParcelType>&, Istream&, bool)"
    );
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::readFields(Cloud<ParcelType>& c)
{
    if (!c.size())
    {
        return;
    }

    Particle<ParcelType>::readFields(c);

    IOField<label> active(c.fieldIOobject("active", IOobject::MUST_READ));
    c.checkFieldIOobject(c, active);

    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeId);

    IOField<scalar>
        nParticle(c.fieldIOobject("nParticle", IOobject::MUST_READ));
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<vector> f(c.fieldIOobject("f", IOobject::MUST_READ));
    c.checkFieldIOobject(c, f);

    IOField<vector> angularMomentum
    (
        c.fieldIOobject("angularMomentum", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, angularMomentum);

    IOField<vector> torque(c.fieldIOobject("torque", IOobject::MUST_READ));
    c.checkFieldIOobject(c, torque);

    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::MUST_READ));
    c.checkFieldIOobject(c, rho);

    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tTurb);

    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::MUST_READ));
    c.checkFieldIOobject(c, UTurb);

    collisionRecordIOList collisionRecords
    (
        IOobject
        (
            "collisionRecords",
            c.time().timeName(),
            c,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    if (collisionRecords.size() != c.size())
    {
        FatalErrorIn
        (
            "void Foam::KinematicParcel<ParcelType>::readFields"
            "(Cloud<ParcelType>& c)"
        )   << "Size of " << collisionRecords.name()
            << " field " << collisionRecords.size()
            << " does not match the number of particles " << c.size()
            << abort(FatalError);
    }

    label i = 0;
    forAllIter(typename Cloud<ParcelType>, c, iter)
    {
        ParcelType& p = iter();

        p.active_ = active[i];
        p.typeId_ = typeId[i];
        p.nParticle_ = nParticle[i];
        p.d_ = d[i];
        p.U_ = U[i];
        p.f_ = f[i];
        p.angularMomentum_ = angularMomentum[i];
        p.rho_ = rho[i];
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];
        p.collisionRecords_ = collisionRecords[i];
        i++;
    }
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::writeFields(const Cloud<ParcelType>& c)
{
    Particle<ParcelType>::writeFields(c);

    label np =  c.size();

    IOField<label> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::NO_READ),
        np
    );
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<vector> f(c.fieldIOobject("f", IOobject::NO_READ), np);
    IOField<vector> angularMomentum
    (
        c.fieldIOobject("angularMomentum", IOobject::NO_READ),
        np
    );
    IOField<vector> torque(c.fieldIOobject("torque", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::NO_READ), np);
    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::NO_READ), np);

    collisionRecordIOList collisionRecords
    (
        IOobject
        (
            "collisionRecords",
            c.time().timeName(),
            c,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        np
    );

    pairDataIOFieldField collisionRecords_pairData
    (
        c.fieldIOobject("collisionRecords_pairData", IOobject::NO_READ),
        np
    );

    label i = 0;

    forAllConstIter(typename Cloud<ParcelType>, c, iter)
    {
        const KinematicParcel<ParcelType>& p = iter();

        active[i] = p.active();
        typeId[i] = p.typeId();
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        U[i] = p.U();
        f[i] = p.f();
        angularMomentum[i] = p.angularMomentum();
        torque[i] = p.torque();
        rho[i] = p.rho();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();
        collisionRecords[i] = p.collisionRecords();
        collisionRecords_pairData[i] = p.collisionRecords().pairData();

        // collisionRecords_pairData[i].setSize
        // (
        //     p.collisionRecords().pairRecords().size()
        // );

        // forAll(p.collisionRecords().pairRecords(), j)
        // {
        //     collisionRecords_pairData[i][j] =
        //         p.collisionRecords().pairRecords()[j].collisionData();
        // }

        i++;
    }

    active.write();
    typeId.write();
    nParticle.write();
    d.write();
    U.write();
    f.write();
    angularMomentum.write();
    torque.write();
    rho.write();
    tTurb.write();
    UTurb.write();
    collisionRecords.write();
    collisionRecords_pairData.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const KinematicParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const Particle<ParcelType>&>(p)
            << token::SPACE << p.active()
            << token::SPACE << p.typeId()
            << token::SPACE << p.nParticle()
            << token::SPACE << p.d()
            << token::SPACE << p.U()
            << token::SPACE << p.f()
            << token::SPACE << p.angularMomentum()
            << token::SPACE << p.torque()
            << token::SPACE << p.rho()
            << token::SPACE << p.tTurb()
            << token::SPACE << p.UTurb()
            << token::SPACE << p.collisionRecords();
    }
    else
    {
        os  << static_cast<const Particle<ParcelType>&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            sizeof(p.active())
          + sizeof(p.typeId())
          + sizeof(p.nParticle())
          + sizeof(p.d())
          + sizeof(p.U())
          + sizeof(p.f())
          + sizeof(p.angularMomentum())
          + sizeof(p.torque())
          + sizeof(p.rho())
          + sizeof(p.tTurb())
          + sizeof(p.UTurb())
        );
        os  << p.collisionRecords();
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const KinematicParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
