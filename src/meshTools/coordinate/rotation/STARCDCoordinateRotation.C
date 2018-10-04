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

#include "STARCDCoordinateRotation.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(STARCDCoordinateRotation, 0);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        STARCDCoordinateRotation,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        STARCDCoordinateRotation,
        objectRegistry
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::vector Foam::STARCDCoordinateRotation::transform(const vector& st) const
{
    return (R_ & st);
}


Foam::vector Foam::STARCDCoordinateRotation::invTransform
(
    const vector& st
) const
{
    return (Rtr_ & st);
}


Foam::tmp<Foam::vectorField> Foam::STARCDCoordinateRotation::transform
(
    const vectorField& st
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::vectorField> Foam::STARCDCoordinateRotation::invTransform
(
    const vectorField& st
) const
{
    NotImplemented;
    return nullptr;
}


const Foam::tensorField& Foam::STARCDCoordinateRotation::Tr() const
{
    NotImplemented;
    return NullObjectRef<tensorField>();
}


Foam::tmp<Foam::tensorField> Foam::STARCDCoordinateRotation::transformTensor
(
    const tensorField& st
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tensor Foam::STARCDCoordinateRotation::transformTensor
(
    const tensor& st
) const
{
    return (R_ & st & Rtr_);
}


Foam::tmp<Foam::tensorField> Foam::STARCDCoordinateRotation::transformTensor
(
    const tensorField& st,
    const labelList& cellMap
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::symmTensorField> Foam::STARCDCoordinateRotation::
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


Foam::symmTensor Foam::STARCDCoordinateRotation::transformVector
(
    const vector& st
) const
{
    return transformPrincipal(R_, st);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tensor Foam::STARCDCoordinateRotation::rotation
(
    const vector& angles,
    const bool degrees
)
{
    scalar z = angles.component(vector::X);    // 1. Rotate about Z
    scalar x = angles.component(vector::Y);    // 2. Rotate about X
    scalar y = angles.component(vector::Z);    // 3. Rotate about Y

    if (degrees)
    {
        x *= degToRad();
        y *= degToRad();
        z *= degToRad();
    }

    return
        tensor
        (
            cos(y)*cos(z) - sin(x)*sin(y)*sin(z),
            -cos(x)*sin(z),
            sin(x)*cos(y)*sin(z) + sin(y)*cos(z),

            cos(y)*sin(z) + sin(x)*sin(y)*cos(z),
            cos(x)*cos(z),
            sin(y)*sin(z) - sin(x)*cos(y)*cos(z),

            -cos(x)*sin(y),
            sin(x),
            cos(x)*cos(y)
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::STARCDCoordinateRotation::STARCDCoordinateRotation()
:
    R_(sphericalTensor::I),
    Rtr_(sphericalTensor::I)
{}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const STARCDCoordinateRotation& r
)
:
    R_(r.R_),
    Rtr_(r.Rtr_)
{}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const vector& rotZrotXrotY,
    const bool degrees
)
:
    R_(rotation(rotZrotXrotY, degrees)),
    Rtr_(R_.T())
{}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const scalar rotZ,
    const scalar rotX,
    const scalar rotY,
    const bool degrees
)
:
    R_(rotation(vector(rotZ, rotX, rotY), degrees)),
    Rtr_(R_.T())
{}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
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


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const dictionary& dict,
    const objectRegistry&
)
:
    STARCDCoordinateRotation(dict)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
