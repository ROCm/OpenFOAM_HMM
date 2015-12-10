/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd
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

#include "solarCalculator.H"
#include "Time.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solarCalculator, 0);

    template<>
    const char* NamedEnum
    <
        solarCalculator::sunDirModel,
        2
    >::names[] =
    {
        "sunDirConstant",
        "sunDirTraking"
    };

    template<>
    const char* NamedEnum
    <
        solarCalculator::sunLModel,
        3
    >::names[] =
    {
        "sunLoadConstant",
        "sunLoadFairWeatherConditions",
        "sunLoadTheoreticalMaximum"
    };
}

const Foam::NamedEnum<Foam::solarCalculator::sunDirModel, 2>
  Foam::solarCalculator::sunDirectionModelTypeNames_;

const Foam::NamedEnum<Foam::solarCalculator::sunLModel, 3>
   Foam::solarCalculator::sunLoadModelTypeNames_;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solarCalculator::calculateBetaTetha()
{
    scalar runTime = 0.0;
    switch (sunDirectionModel_)
    {
        case mSunDirTraking:
        {
            runTime = mesh_.time().value();
            break;
        }
        case mSunDirConstant:
        {
            break;
        }
    }

    scalar LSM = 15.0*(readScalar(dict_.lookup("localStandardMeridian")));

    scalar D = readScalar(dict_.lookup("startDay")) + runTime/86400.0;
    scalar M = 6.24004 + 0.0172*D;
    scalar EOT = -7.659*sin(M) + 9.863*sin(2*M + 3.5932);

    startTime_ = readScalar(dict_.lookup("startTime"));
    scalar LST =  startTime_ + runTime/3600.0;

    scalar LON = readScalar(dict_.lookup("longitude"));

    scalar AST = LST + EOT/60.0 + (LON - LSM)/15;

    scalar delta = 23.45*sin(degToRad((360*(284 + D))/365));

    scalar H = degToRad(15*(AST - 12));

    scalar L = degToRad(readScalar(dict_.lookup("latitude")));

    scalar deltaRad = degToRad(delta);
    beta_ = max(asin(cos(L)*cos(deltaRad)*cos(H) + sin(L)*sin(deltaRad)), 1e-3);
    tetha_ = acos((sin(beta_)*sin(L) - sin(deltaRad))/(cos(beta_)*cos(L)));

    if (debug)
    {
        Info << tab << "altitude : " << radToDeg(beta_) << endl;
        Info << tab << "azimuth  : " << radToDeg(tetha_) << endl;
    }
}


void Foam::solarCalculator::calculateSunDirection()
{

    dict_.lookup("gridUp") >> gridUp_;
    gridUp_ /= mag(gridUp_);

    dict_.lookup("gridEast") >> eastDir_;
    eastDir_ /= mag(eastDir_);

    coord_.reset
    (
        new coordinateSystem("grid", vector::zero, gridUp_, eastDir_)
    );

    direction_.z() = -sin(beta_);
    direction_.y() =  cos(beta_)*cos(tetha_); //North
    direction_.x() =  cos(beta_)*sin(tetha_); //East

    direction_ /= mag(direction_);

    if (debug)
    {
        Info<< "Sun direction in absolute coordinates : " << direction_ <<endl;
    }

    direction_ = coord_->R().transform(direction_);

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
            if (dict_.found("sunDirection"))
            {
                dict_.lookup("sunDirection") >> direction_;
                direction_ /= mag(direction_);
            }
            else
            {
                calculateBetaTetha();
                calculateSunDirection();
            }

            break;
        }
        case mSunDirTraking:
        {
            if (word(mesh_.ddtScheme("default")) == "steadyState")
            {
                FatalErrorInFunction
                    << " Sun direction model can not be sunDirtracking if the "
                    << " case is steady " << nl << exit(FatalError);
            }

            dict_.lookup("sunTrackingUpdateInterval") >>
                sunTrackingUpdateInterval_;

            calculateBetaTetha();
            calculateSunDirection();
            break;
        }
    }

    switch (sunLoadModel_)
    {
        case mSunLoadConstant:
        {
            dict_.lookup("directSolarRad") >> directSolarRad_;
            dict_.lookup("diffuseSolarRad") >> diffuseSolarRad_;
            break;
        }
        case mSunLoadFairWeatherConditions:
        {
            A_ = readScalar(dict_.lookup("A"));
            B_ = readScalar(dict_.lookup("B"));

            if (dict_.found("beta"))
            {
                dict_.lookup("beta") >> beta_;
            }
            else
            {
                calculateBetaTetha();
            }

            directSolarRad_ = A_/exp(B_/sin(beta_));

            groundReflectivity_ =
                readScalar(dict_.lookup("groundReflectivity"));

            break;
        }
        case mSunLoadTheoreticalMaximum:
        {
            Setrn_ = readScalar(dict_.lookup("Setrn"));
            SunPrime_ = readScalar(dict_.lookup("SunPrime"));
            directSolarRad_ = Setrn_*SunPrime_;

            groundReflectivity_ =
                readScalar(dict_.lookup("groundReflectivity"));
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
    direction_(vector::zero),
    directSolarRad_(0.0),
    diffuseSolarRad_(0.0),
    groundReflectivity_(0.0),
    A_(0.0),
    B_(0.0),
    beta_(0.0),
    tetha_(0.0),
    Setrn_(0.0),
    SunPrime_(0.0),
    C_(readScalar(dict.lookup("C"))),
    sunDirectionModel_
    (
        sunDirectionModelTypeNames_.read(dict.lookup("sunDirectionModel"))
    ),
    sunLoadModel_
    (
        sunLoadModelTypeNames_.read(dict.lookup("sunLoadModel"))
    ),
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
        case mSunDirTraking:
        {
            calculateBetaTetha();
            calculateSunDirection();
            directSolarRad_ = A_/exp(B_/sin(max(beta_, ROOTVSMALL)));
            break;
        }
    }
}


// ************************************************************************* //
