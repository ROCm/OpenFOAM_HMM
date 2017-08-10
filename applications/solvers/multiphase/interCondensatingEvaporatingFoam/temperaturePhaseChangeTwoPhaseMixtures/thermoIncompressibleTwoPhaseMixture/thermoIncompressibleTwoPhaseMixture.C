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

#include "thermoIncompressibleTwoPhaseMixture.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoIncompressibleTwoPhaseMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermoIncompressibleTwoPhaseMixture::thermoIncompressibleTwoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    incompressibleTwoPhaseMixture(U, phi),

    kappa1_
    (
        "kappa1",
        dimEnergy/dimTime/dimLength/dimTemperature,
        subDict(phase1Name_).lookup("kappa")
    ),
    kappa2_
    (
        "kappa2",
        kappa1_.dimensions(),
        subDict(phase2Name_).lookup("kappa")
    ),

    Cp1_
    (
        "Cp1",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase1Name_).lookup("Cp")
    ),
    Cp2_
    (
        "Cp2",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase2Name_).lookup("Cp")
    ),

    Cv1_
    (
        "Cv1",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase1Name_).lookup("Cv")
    ),

    Cv2_
    (
        "Cv2",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase2Name_).lookup("Cv")
    ),

    Hf1_
    (
        "Hf1",
        dimEnergy/dimMass,
        subDict(phase1Name_).lookup("hf")
    ),

    Hf2_
    (
        "Hf2",
        dimEnergy/dimMass,
        subDict(phase2Name_).lookup("hf")
    )
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool Foam::thermoIncompressibleTwoPhaseMixture::read()
{
    if (incompressibleTwoPhaseMixture::read())
    {
        subDict(phase1Name_).lookup("kappa") >> kappa1_;
        subDict(phase2Name_).lookup("kappa") >> kappa2_;

        subDict(phase1Name_).lookup("Cp") >> Cp1_;
        subDict(phase2Name_).lookup("Cp") >> Cp2_;

        subDict(phase1Name_).lookup("Cv") >> Cv1_;
        subDict(phase2Name_).lookup("Cv") >> Cv2_;

        subDict(phase1Name_).lookup("hf") >> Hf1_;
        subDict(phase2Name_).lookup("hf") >> Hf2_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
