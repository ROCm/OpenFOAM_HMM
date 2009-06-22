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

#include "Particle.H"
#include "IOstreams.H"
#include "IOPosition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
Foam::string Foam::Particle<ParticleType>::propHeader =
    "(Px Py Pz) cellI origProc origId";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Particle<ParticleType>::Particle
(
    const Cloud<ParticleType>& cloud,
    Istream& is,
    bool readFields
)
:
    cloud_(cloud),
    facei_(-1),
    stepFraction_(0.0),
    origProc_(Pstream::myProcNo()),
    origId_(-1)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> position_ >> celli_ >> origProc_ >> origId_;
    }
    else
    {
        // In binary read all particle data - needed for parallel transfer
        is.read
        (
            reinterpret_cast<char*>(&position_),
            sizeof(position_)
            + sizeof(celli_)
            + sizeof(facei_)
            + sizeof(stepFraction_)
            + sizeof(origProc_)
            + sizeof(origId_)
        );
    }

    if (celli_ == -1)
    {
        celli_ = cloud_.pMesh().findCell(position_);
    }

    // Check state of Istream
    is.check("Particle<ParticleType>::Particle(Istream&)");
}


template<class ParticleType>
void Foam::Particle<ParticleType>::writeFields
(
    const Cloud<ParticleType>& c
)
{
    // Write the cloud position file
    IOPosition<ParticleType> ioP(c);
    ioP.write();
}


template<class ParticleType>
Foam::Ostream& Foam::operator<<(Ostream& os, const Particle<ParticleType>& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os << p.position_
           << token::SPACE << p.celli_
           << token::SPACE << p.origProc_
           << token::SPACE << p.origId_;
    }
    else
    {
        // In binary write both celli_ and facei_, needed for parallel transfer
        os.write
        (
            reinterpret_cast<const char*>(&p.position_),
            sizeof(p.position_)
          + sizeof(p.celli_)
          + sizeof(p.facei_)
          + sizeof(p.stepFraction_)
          + sizeof(p.origProc_)
          + sizeof(p.origId_)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const Particle<ParticleType>&)");

    return os;
}


// ************************************************************************* //
