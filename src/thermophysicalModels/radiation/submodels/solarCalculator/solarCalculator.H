/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenCFD Ltd.
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

Class
    Foam::solarCalculator

Description
    The solar calculator model provides information about the Sun direction
    and Sun load model. The available models are:

    For the Sun direction:
    1) SunDirConstant : the direction is given in 'sunDirection'
    2) SunDirTracking : the direction is calculated from the following
       parameters:
            localStandardMeridian : GMT (Local Zone Meridian) in hours
            startDay :  day from 1 to 365)
            startTime:  in hours
            longitude:  in degrees
            latitude:   in degrees
            gridUp:     grid orientation upwards
            gridEast    grid orientation eastwards

        This model should be use in transient calculations.
        The keyword 'sunTrackingUpdateInterval' (in hours) specifies on which
        interval is the Sun direction updated.


    Solar Load models available:
    1) SunLoadConstant: direct and diffusive heat fluxes are provided by the
       entries 'directSolarRad' and 'diffuseSolarRad'

    2) SunLoadFairWeatherConditions: The solar fluxes are calculated following
       the Fair Weather Conditions Method from the ASHRAE Handbook. The entries
       are:
            skyCloudCoverFraction: Fraction of sky covered by clouds (0-1)
            A :     Apparent solar irradiation at air mass m = 0
            B :     Atmospheric extinction coefficient
            beta:   Solar altitude (in degrees) above the horizontal. This
            can be read or calculated providing the respective parameters
            for Sun position explained above.
            groundReflectivity : ground reflectivity

            In this model the flux is calculated as:

                directSolarRad =
                    (1 - 0.75*skyCloudCoverFraction^3)*A/exp(B/sin(beta));

    3) SunLoadTheoreticalMaximum: The entries are:
            Setrn
            SunPrime:
            groundReflectivity : ground reflectivity

            In this model the flux is calculated as:

                directSolarRad = Setrn*SunPrime;

    The diffuse on vertical/horizontal walls and ground-reflected radiation are
    calculated following the ASHRAE Handbook.


SourceFiles
    solarCalculator.C

\*---------------------------------------------------------------------------*/

#ifndef solarCalculator_H
#define solarCalculator_H

#include "fvMesh.H"
#include "meshTools.H"
#include "DynamicField.H"
#include "HashSet.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solarCalculator Declaration
\*---------------------------------------------------------------------------*/

class solarCalculator
{
public:

    //  Public enumeration

        //- Sun direction models
        enum sunDirModel
        {
            mSunDirConstant,
            mSunDirTracking
        };

        //- Direct sun load models
        enum sunLModel
        {
            mSunLoadConstant,
            mSunLoadFairWeatherConditions,
            mSunLoadTheoreticalMaximum
        };


protected:

         //- Sun direction models
        static const Enum<sunDirModel> sunDirectionModelTypeNames_;

        //- Sun load models
        static const Enum<sunLModel> sunLoadModelTypeNames_;

private:

    // Private data


        //- Reference to mesh
        const fvMesh& mesh_;

        //- Dictionary
        dictionary dict_;

        //- Direction
        vector direction_;

        //- Direct solar irradiation
        scalar directSolarRad_;

        //- Diffuse solar irradiation on vertical surfaces
        scalar diffuseSolarRad_;

        //- Ground reflectivity
        scalar groundReflectivity_;

        //- Fair weather direct solar load model parameters
        scalar A_;
        scalar B_;
        scalar beta_;
        scalar tetha_;

        //- Sky cloud cover fraction [0-1]
        scalar skyCloudCoverFraction_;


        //- Maximum theoretical direct solar load model parameters
        scalar Setrn_;
        scalar SunPrime_;


        //- Diffusive solar load model parameters
        scalar C_;

        //- Sun direction model
        sunDirModel sunDirectionModel_;

        //- Sun load model
        sunLModel sunLoadModel_;

        //- Grid coordinate system
        autoPtr<coordinateSystem> coord_;

        //- East grid orientation
        vector eastDir_;

        //- Up grid orientation
        vector gridUp_;

        //- Interval in decimal hours to update Sun direction for SunDirTracking
        scalar sunTrackingUpdateInterval_;

        //- Start time for the Sun position (decimal hours)
        scalar startTime_;


        //- No copy construct
        solarCalculator(const solarCalculator&) = delete;

        //- No copy assignment
        void operator=(const solarCalculator&) = delete;


    // Private members

        //- Init
        void init();

        //- Calculate beta and tetha angles
        void calculateBetaTetha();

        //- Calculate Sun direction
        void calculateSunDirection();


public:

    // Declare name of the class and its debug switch
    ClassName("solarCalculator");


    // Constructors

        //- Construct from dictionary
        solarCalculator(const dictionary&, const fvMesh&);


    //- Destructor
    ~solarCalculator();


    // Member Functions

        // Access

            //- const acess to direction
            const vector direction() const
            {
                return direction_;
            }

            //- Non-const access to direction
            vector& direction()
            {
                return direction_;
            }

             //- Return direct solar irradiation
            scalar directSolarRad()
            {
                return directSolarRad_;
            }

            //- Return diffuse solar irradiation
            scalar diffuseSolarRad()
            {
                return diffuseSolarRad_;
            }

            //- Return C constant
            scalar C()
            {
                return C_;
            }

            //- Return beta
            scalar beta()
            {
                return beta_;
            }

            //- Return tetha
            scalar tetha()
            {
                return tetha_;
            }

            //- Return Sun direction model
            sunDirModel sunDirectionModel() const
            {
                return sunDirectionModel_;
            }

            //- Return Sun load model
            sunLModel sunLoadModel() const
            {
                return sunLoadModel_;
            }

            //- Return groundReflectivity
            scalar groundReflectivity()
            {
                return groundReflectivity_;
            }

            //- Return coordinateSystem
            const coordinateSystem& coord()
            {
                return *coord_;
            }

            //- Return sunTrackingUpdateInterval
            scalar sunTrackingUpdateInterval()
            {
                return sunTrackingUpdateInterval_;
            }

            //- Return startTime
            scalar startTime()
            {
                return startTime_;
            }


        //- Recalculate
        void correctSunDirection();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
