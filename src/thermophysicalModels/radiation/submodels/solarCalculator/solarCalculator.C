/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd
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

#include "solarCalculator.H"
#include "Time.H"
#include "unitConversion.H"
#include "constants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solarCalculator, 0);
}


const Foam::Enum
<
    Foam::solarCalculator::sunDirModel
>
Foam::solarCalculator::sunDirectionModelTypeNames_
({
    { sunDirModel::mSunDirConstant, "sunDirConstant" },
    { sunDirModel::mSunDirTracking, "sunDirTracking" },
});


const Foam::Enum
<
    Foam::solarCalculator::sunLModel
>
Foam::solarCalculator::sunLoadModelTypeNames_
({
    { sunLModel::mSunLoadConstant, "sunLoadConstant" },
    {
        sunLModel::mSunLoadFairWeatherConditions,
        "sunLoadFairWeatherConditions"
    },
    { sunLModel::mSunLoadTheoreticalMaximum, "sunLoadTheoreticalMaximum" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solarCalculator::calculateBetaTetha()
{
    scalar runTime = 0.0;
    switch (sunDirectionModel_)
    {
        case mSunDirTracking:
        {
            runTime = mesh_.time().value();
            break;
        }
        case mSunDirConstant:
        {
            break;
        }
    }

    scalar LSM = 15.0*(dict_.get<scalar>("localStandardMeridian"));

    scalar D = dict_.get<scalar>("startDay") + runTime/86400.0;
    scalar M = 6.24004 + 0.0172*D;
    scalar EOT = -7.659*sin(M) + 9.863*sin(2*M + 3.5932);

    dict_.readEntry("startTime", startTime_);

    scalar LST =  startTime_ + runTime/3600.0;

    scalar LON = dict_.get<scalar>("longitude");

    scalar AST = LST + EOT/60.0 + (LON - LSM)/15;

    scalar delta = 23.45*sin(degToRad((360*(284 + D))/365));

    scalar H = degToRad(15*(AST - 12));

    scalar L = degToRad(dict_.get<scalar>("latitude"));

    scalar deltaRad = degToRad(delta);
    beta_ = max(asin(cos(L)*cos(deltaRad)*cos(H) + sin(L)*sin(deltaRad)), 1e-3);
    tetha_ = acos((sin(beta_)*sin(L) - sin(deltaRad))/(cos(beta_)*cos(L)));

    // theta is the angle between the SOUTH axis and the Sun
    // If the hour angle is lower than zero (morning) the Sun is positioned
    // on the East side.
    if (H < 0)
    {
        tetha_ += 2*(constant::mathematical::pi - tetha_);
    }

    if (debug)
    {
        Info << tab << "altitude : " << radToDeg(beta_) << endl;
        Info << tab << "azimuth  : " << radToDeg(tetha_) << endl;
    }
}


void Foam::solarCalculator::calculateSunDirection()
{
    gridUp_  = normalised(dict_.get<vector>("gridUp"));
    eastDir_ = normalised(dict_.get<vector>("gridEast"));

    coord_.reset
    (
        new coordinateSystem("grid", Zero, gridUp_, eastDir_)
    );

    // Assuming 'z' vertical, 'y' North and 'x' East
    direction_.z() = -sin(beta_);
    direction_.y() =  cos(beta_)*cos(tetha_); // South axis
    direction_.x() =  cos(beta_)*sin(tetha_); // West axis

    direction_.normalise();

    if (debug)
    {
        Info<< "Sun direction in absolute coordinates : " << direction_ <<endl;
    }

    // Transform to actual coordinate system
    direction_ = coord_->transform(direction_);

    if (debug)
    {
        Info<< "Sun direction in the Grid coordinates : " << direction_ <<endl;
    }
}


void Foam::solarCalculator::init()
{
    switch (sunDirectionModel_)
    {
        case mSunDirConstant:
        {
            if (dict_.readIfPresent("sunDirection", direction_))
            {
                direction_.normalise();
            }
            else
            {
                calculateBetaTetha();
                calculateSunDirection();
            }

            break;
        }
        case mSunDirTracking:
        {
            if (word(mesh_.ddtScheme("default")) == "steadyState")
            {
                FatalErrorInFunction
                    << " Sun direction model can not be sunDirtracking if the "
                    << " case is steady " << nl << exit(FatalError);
            }

            dict_.readEntry
            (
                "sunTrackingUpdateInterval",
                sunTrackingUpdateInterval_
            );

            calculateBetaTetha();
            calculateSunDirection();
            break;
        }
    }

    switch (sunLoadModel_)
    {
        case mSunLoadConstant:
        {
            dict_.readEntry("directSolarRad", directSolarRad_);
            dict_.readEntry("diffuseSolarRad", diffuseSolarRad_);
            break;
        }
        case mSunLoadFairWeatherConditions:
        {
            dict_.readIfPresent
            (
                "skyCloudCoverFraction",
                skyCloudCoverFraction_
            );

            dict_.readEntry("A", A_);
            dict_.readEntry("B", B_);

            if (!dict_.readIfPresent("beta", beta_))
            {
                calculateBetaTetha();
            }

            directSolarRad_ =
                (1.0 - 0.75*pow(skyCloudCoverFraction_, 3.0))
              * A_/exp(B_/sin(beta_));

            dict_.readEntry("groundReflectivity", groundReflectivity_);
            break;
        }
        case mSunLoadTheoreticalMaximum:
        {
            dict_.readEntry("Setrn", Setrn_);
            dict_.readEntry("SunPrime", SunPrime_);
            directSolarRad_ = Setrn_*SunPrime_;

            dict_.readEntry("groundReflectivity", groundReflectivity_);
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solarCalculator::solarCalculator
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_(dict),
    direction_(Zero),
    directSolarRad_(0.0),
    diffuseSolarRad_(0.0),
    groundReflectivity_(0.0),
    A_(0.0),
    B_(0.0),
    beta_(0.0),
    tetha_(0.0),
    skyCloudCoverFraction_(0.0),
    Setrn_(0.0),
    SunPrime_(0.0),
    C_(dict.get<scalar>("C")),
    sunDirectionModel_
    (
        sunDirectionModelTypeNames_.get("sunDirectionModel", dict)
    ),
    sunLoadModel_(sunLoadModelTypeNames_.get("sunLoadModel", dict)),
    coord_()
{
    init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solarCalculator::~solarCalculator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solarCalculator::correctSunDirection()
{
    switch (sunDirectionModel_)
    {
        case mSunDirConstant:
        {
            break;
        }
        case mSunDirTracking:
        {
            calculateBetaTetha();
            calculateSunDirection();
            directSolarRad_ = A_/exp(B_/sin(max(beta_, ROOTVSMALL)));
            break;
        }
    }
}


// ************************************************************************* //
