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
    namespace coordinateRotations
    {
        defineTypeName(euler);

        // Standard short name
        addNamedToRunTimeSelectionTable
        (
            coordinateRotation,
            euler,
            dictionary,
            euler
        );

        // Longer name - Compat 1806
        addNamedToRunTimeSelectionTable
        (
            coordinateRotation,
            euler,
            dictionary,
            EulerRotation
        );
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tensor Foam::coordinateRotations::euler::rotation
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateRotations::euler::euler()
:
    coordinateRotation(),
    angles_(Zero),
    degrees_(true)
{}


Foam::coordinateRotations::euler::euler(const euler& crot)
:
    coordinateRotation(crot),
    angles_(crot.angles_),
    degrees_(crot.degrees_)
{}


Foam::coordinateRotations::euler::euler
(
    const vector& phiThetaPsi,
    bool degrees
)
:
    coordinateRotation(),
    angles_(phiThetaPsi),
    degrees_(degrees)
{}


Foam::coordinateRotations::euler::euler
(
    scalar phi,
    scalar theta,
    scalar psi,
    bool degrees
)
:
    coordinateRotation(),
    angles_(phi, theta, psi),
    degrees_(degrees)
{}


Foam::coordinateRotations::euler::euler(const dictionary& dict)
:
    coordinateRotation(),
    angles_(dict.getCompat<vector>("angles", {{"rotation", 1806}})),
    degrees_(dict.lookupOrDefault("degrees", true))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::coordinateRotations::euler::clear()
{
    angles_ = Zero;
    degrees_ = true;
}


Foam::tensor Foam::coordinateRotations::euler::R() const
{
    return euler::rotation(angles_, degrees_);
}


void Foam::coordinateRotations::euler::write(Ostream& os) const
{
    os  << "euler-angles(" << (degrees_ ? "deg" : "rad") << "): " << angles_;
}


void Foam::coordinateRotations::euler::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.beginBlock(keyword);

    os.writeEntry("type", type());
    os.writeEntry("angles", angles_);
    if (!degrees_)
    {
        os.writeEntry("degrees", "false");
    }

    os.endBlock();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
