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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::tensor
Foam::sixDofRigidBodyMotion::rotationTensorX(scalar phi) const
{
    return tensor
    (
        1, 0, 0,
        0, Foam::cos(phi), -Foam::sin(phi),
        0, Foam::sin(phi), Foam::cos(phi)
    );
}


inline Foam::tensor
Foam::sixDofRigidBodyMotion::rotationTensorY(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), 0, Foam::sin(phi),
        0, 1, 0,
        -Foam::sin(phi), 0, Foam::cos(phi)
    );
}


inline Foam::tensor
Foam::sixDofRigidBodyMotion::rotationTensorZ(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), -Foam::sin(phi), 0,
        Foam::sin(phi), Foam::cos(phi), 0,
        0, 0, 1
    );
}


inline const Foam::point& Foam::sixDofRigidBodyMotion::centreOfMass() const
{
    return centreOfMass_;
}


inline const Foam::point& Foam::sixDofRigidBodyMotion::refCentreOfMass() const
{
    return refCentreOfMass_;
}


inline const Foam::diagTensor&
Foam::sixDofRigidBodyMotion::momentOfInertia() const
{
    return momentOfInertia_;
}


inline Foam::scalar Foam::sixDofRigidBodyMotion::mass() const
{
    return mass_;
}


inline const Foam::tensor& Foam::sixDofRigidBodyMotion::Q() const
{
    return Q_;
}


inline const Foam::vector& Foam::sixDofRigidBodyMotion::v() const
{
    return v_;
}


inline const Foam::vector& Foam::sixDofRigidBodyMotion::a() const
{
    return a_;
}


inline const Foam::vector& Foam::sixDofRigidBodyMotion::pi() const
{
    return pi_;
}


inline const Foam::vector& Foam::sixDofRigidBodyMotion::tau() const
{
    return tau_;
}


inline Foam::point& Foam::sixDofRigidBodyMotion::centreOfMass()
{
    return centreOfMass_;
}


inline Foam::point& Foam::sixDofRigidBodyMotion::refCentreOfMass()
{
    return refCentreOfMass_;
}


inline Foam::diagTensor& Foam::sixDofRigidBodyMotion::momentOfInertia()
{
    return momentOfInertia_;
}


inline Foam::scalar& Foam::sixDofRigidBodyMotion::mass()
{
    return mass_;
}


inline Foam::tensor& Foam::sixDofRigidBodyMotion::Q()
{
    return Q_;
}


inline Foam::vector& Foam::sixDofRigidBodyMotion::v()
{
    return v_;
}


inline Foam::vector& Foam::sixDofRigidBodyMotion::a()
{
    return a_;
}


inline Foam::vector& Foam::sixDofRigidBodyMotion::pi()
{
    return pi_;
}


inline Foam::vector& Foam::sixDofRigidBodyMotion::tau()
{
    return tau_;
}


// ************************************************************************* //
