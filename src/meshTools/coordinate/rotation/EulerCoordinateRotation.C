/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "EulerCoordinateRotation.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(EulerCoordinateRotation, 0);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        EulerCoordinateRotation,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        EulerCoordinateRotation,
        objectRegistry
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tensor Foam::EulerCoordinateRotation::rotation
(
    const vector& angles,
    bool degrees
)
{
    scalar phi   = angles.component(vector::X); // 1. Rotate about Z
    scalar theta = angles.component(vector::Y); // 2. Rotate about X
    scalar psi   = angles.component(vector::Z); // 3. Rotate about Z

    if (degrees)
    {
        phi   *= degToRad();
        theta *= degToRad();
        psi   *= degToRad();
    }

    const scalar c1 = cos(phi);   const scalar s1 = sin(phi);
    const scalar c2 = cos(theta); const scalar s2 = sin(theta);
    const scalar c3 = cos(psi);   const scalar s3 = sin(psi);

    // Compare
    // https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
    //
    // Z1-X2-Z3 rotation

    return
        tensor
        (
            c1*c3 - c2*s1*s3, -c1*s3 - c2*c3*s1,  s1*s2,
            c3*s1 + c1*c2*s3,  c1*c2*c3 - s1*s3, -c1*s2,
            s2*s3,             c3*s2,             c2
        );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::vector Foam::EulerCoordinateRotation::transform(const vector& st) const
{
    return (R_ & st);
}


Foam::vector Foam::EulerCoordinateRotation::invTransform
(
    const vector& st
) const
{
    return (Rtr_ & st);
}


Foam::tmp<Foam::vectorField> Foam::EulerCoordinateRotation::transform
(
    const vectorField& st
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::vectorField> Foam::EulerCoordinateRotation::invTransform
(
    const vectorField& st
) const
{
    NotImplemented;
    return nullptr;
}


const Foam::tensorField& Foam::EulerCoordinateRotation::Tr() const
{
    NotImplemented;
    return NullObjectRef<tensorField>();
}


Foam::tmp<Foam::tensorField> Foam::EulerCoordinateRotation::transformTensor
(
    const tensorField& st
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tensor Foam::EulerCoordinateRotation::transformTensor
(
    const tensor& st
) const
{
    return (R_ & st & Rtr_);
}


Foam::tmp<Foam::tensorField> Foam::EulerCoordinateRotation::transformTensor
(
    const tensorField& st,
    const labelList& cellMap
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::symmTensorField> Foam::EulerCoordinateRotation::
transformVector
(
    const vectorField& st
) const
{
    tmp<symmTensorField> tfld(new symmTensorField(st.size()));
    symmTensorField& fld = tfld.ref();

    forAll(fld, i)
    {
        fld[i] = transformPrincipal(R_, st[i]);
    }
    return tfld;
}


Foam::symmTensor Foam::EulerCoordinateRotation::transformVector
(
    const vector& st
) const
{
    return transformPrincipal(R_, st);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EulerCoordinateRotation::EulerCoordinateRotation()
:
    R_(sphericalTensor::I),
    Rtr_(sphericalTensor::I)
{}


Foam::EulerCoordinateRotation::EulerCoordinateRotation
(
    const EulerCoordinateRotation& r
)
:
    R_(r.R_),
    Rtr_(r.Rtr_)
{}


Foam::EulerCoordinateRotation::EulerCoordinateRotation
(
    const vector& phiThetaPsi,
    const bool degrees
)
:
    R_(rotation(phiThetaPsi, degrees)),
    Rtr_(R_.T())
{}


Foam::EulerCoordinateRotation::EulerCoordinateRotation
(
    const scalar phi,
    const scalar theta,
    const scalar psi,
    const bool degrees
)
:
    R_(rotation(vector(phi, theta, psi), degrees)),
    Rtr_(R_.T())
{}


Foam::EulerCoordinateRotation::EulerCoordinateRotation
(
    const dictionary& dict
)
:
    R_
    (
        rotation
        (
            dict.get<vector>("rotation"),
            dict.lookupOrDefault("degrees", true)
        )
    ),
    Rtr_(R_.T())
{}


Foam::EulerCoordinateRotation::EulerCoordinateRotation
(
    const dictionary& dict,
    const objectRegistry&
)
:
    EulerCoordinateRotation(dict)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
