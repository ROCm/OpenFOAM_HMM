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

#include "ParticleTrackingData.H"

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::ParticleTrackingData<ParcelType>::readProperties
(
    const Cloud<ParcelType>& cloud
)
{
    IOobject propsDictHeader
    (
        "particleTrackingProperties",
        cloud.db().time().timeName(),
        "uniform/Lagrangian"/cloud.name(),
        cloud.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (propsDictHeader.headerOk())
    {
        const IOdictionary propsDict(propsDictHeader);

        word procName("processor" + name(Pstream::myProcNo()));
        if (propsDict.found(procName))
        {
            propsDict.subDict(procName).lookup("particleCount") >>
                PARTICLE_COUNT;
        }
    }
}


template<class ParcelType>
void Foam::ParticleTrackingData<ParcelType>::writeProperties
(
    const Cloud<ParcelType>& cloud
)
{
    if (cloud.db().time().outputTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                "particleTrackingProperties",
                cloud.db().time().timeName(),
                "uniform/Lagrangian"/cloud.name(),
                cloud.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        word procName("processor" + name(Pstream::myProcNo()));
        propsDict.add(procName, dictionary());
        propsDict.subDict(procName).add("particleCount", PARTICLE_COUNT);

        propsDict.regIOobject::write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ParticleTrackingData<ParcelType>::ParticleTrackingData
(
    const Cloud<ParcelType>& cloud,
    Istream& is,
    bool readFields
)
:
    cloud_(cloud),
    origProc_(-1),
    id_(-1)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> origProc_ >> id_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&origProc_),
                sizeof(origProc_) + sizeof(id_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "ParticleTrackingData<ParcelType>::ParticleTrackingData"
        "("
            "Istream&, "
            "bool"
        ")"
    );
}


template<class ParcelType>
void Foam::ParticleTrackingData<ParcelType>::readFields
(
    Cloud<ParcelType>& c
)
{
    if (!c.size())
    {
        return;
    }

    readProperties(c);

    IOField<label> origProc(c.fieldIOobject("origProc", IOobject::MUST_READ));
    c.checkFieldIOobject(c, origProc);

    IOField<label> id(c.fieldIOobject("id", IOobject::MUST_READ));
    c.checkFieldIOobject(c, id);

    label i = 0;
    forAllIter(typename Cloud<ParcelType>, c, iter)
    {
        ParcelType& p = iter();
        p.origProc_ = origProc[i];
        p.id_ = id[i];
        i++;
    }
}


template<class ParcelType>
void Foam::ParticleTrackingData<ParcelType>::writeFields
(
    const Cloud<ParcelType>& c
)
{
    writeProperties(c);

    const label np = c.size();

    IOField<label> origProc
        (
            c.fieldIOobject("origProc", IOobject::NO_READ),
            np
        );
    IOField<label> id(c.fieldIOobject("id", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<ParcelType>, c, iter)
    {
        const ParcelType& p = iter();

        origProc[i] = p.origProc();
        id[i] = p.id();
        i++;
    }

    origProc.write();
    id.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ParticleTrackingData<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << p.origProc_ << token::SPACE << p.id_ << token::SPACE;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&p.origProc_),
            sizeof(p.origProc_) + sizeof(p.id_)
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<"
        "("
            "Ostream&, "
            "const ParticleTrackingData<ParcelType>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
