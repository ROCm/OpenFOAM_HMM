/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "sixDofRigidBodyMotionState.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDofRigidBodyMotionState::write(Ostream& os) const
{
    os.writeKeyword("centreOfMass")
        << centreOfMass_ << token::END_STATEMENT << nl;
    os.writeKeyword("Q")
        << Q_ << token::END_STATEMENT << nl;
    os.writeKeyword("v")
        << v_ << token::END_STATEMENT << nl;
    os.writeKeyword("a")
        << a_ << token::END_STATEMENT << nl;
    os.writeKeyword("pi")
        << pi_ << token::END_STATEMENT << nl;
    os.writeKeyword("tau")
        << tau_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is, sixDofRigidBodyMotionState& sDofRBMS
)
{
    is  >> sDofRBMS.centreOfMass_
        >> sDofRBMS.Q_
        >> sDofRBMS.v_
        >> sDofRBMS.a_
        >> sDofRBMS.pi_
        >> sDofRBMS.tau_;

    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::sixDofRigidBodyMotionState&)"
    );

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sixDofRigidBodyMotionState& sDofRBMS
)
{
    os  << token::SPACE << sDofRBMS.centreOfMass()
        << token::SPACE << sDofRBMS.Q()
        << token::SPACE << sDofRBMS.v()
        << token::SPACE << sDofRBMS.a()
        << token::SPACE << sDofRBMS.pi()
        << token::SPACE << sDofRBMS.tau();

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::sixDofRigidBodyMotionState&)"
    );

    return os;
}


// ************************************************************************* //
