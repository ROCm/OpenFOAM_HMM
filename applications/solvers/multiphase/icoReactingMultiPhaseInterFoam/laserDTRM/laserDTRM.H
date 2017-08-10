/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd
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
    Foam::radiation::laserDTRM

Description
    Works well for combustion applications where optical thickness, tau is
    large, i.e. tau = a*L > 3 (L = distance between objects)

    Assumes
     - all surfaces are diffuse
     - tends to over predict radiative fluxes from sources/sinks
       *** SOURCES NOT CURRENTLY INCLUDED ***

SourceFiles
    laserDTRM.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModellaserDTRM_H
#define radiationModellaserDTRM_H

#include "radiationModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Cloud.H"
#include "DTRMParticle.H"
#include "reflectionModel.H"
#include "Function1.H"
#include "interpolation2DTable.H"
#include "labelField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class laserDTRM Declaration
\*---------------------------------------------------------------------------*/

class laserDTRM
:
    public radiationModel
{
public:

    //  Public enumeration

        //- Modes of power distribution
        enum powerDistributionMode
        {
            pdGaussian,
            pdManual,
            pdUniform
        };

private:


    // Private data


        static const NamedEnum<powerDistributionMode, 3> powerDistypeNames_;

        //- Operating mode for power distribution
        powerDistributionMode mode_;

        //- DTRM particle cloud
        Cloud<DTRMParticle> DTRMCloud_;

        //- Number of particles
        label nParticles_;

        //- Number of theta angles
        label ndTheta_;

        //- Number of radials
        label ndr_;

        //- Maximum tracking length for particles
        scalar maxTrackLength_;

        // Laser parameters

            //- Focal laser position
            autoPtr<Function1<point>> focalLaserPosition_;

            //- Direction
            autoPtr<Function1<vector>> laserDirection_;

            //- Focal lase radius
            scalar focalLaserRadius_;

            //- Quality beam laser (tan of the angle between the main direction
            //  (laserDirection) and the spreading rays (0: flat, 1:fully)
            scalar qualityBeamLaser_;


            // Gaussian power distribution entries

                //- Sigma deviation
                scalar sigma_;

                //- Power
                autoPtr<Function1<scalar>> laserPower_;


            // Manual power distribution entries

                //- Manual power distribution table. (theta, radius)
                autoPtr<interpolation2DTable<scalar>> powerDistribution_;


            // Reflection sub-model

                //- Reflection model
                autoPtr<reflectionModel> reflection_;

                //- Reflection switch
                bool reflectionSwitch_;

                //- Phase in which the particles are inserted
                word initialPhase_;

                //- Phase name for absorbing medium 1
                word alpha1_;

                //- Phase name for absorbing medium 2
                word alpha2_;


        // Fields


            //- Incident radiative heat flux [W/m2]
            volScalarField Qin_;

            //- Absorption coefficient
            volScalarField a_;

            //- Emission coefficient
            volScalarField e_;

            //- Emission contribution
            volScalarField E_;

            //- Volumetric heat source [W/m3]
            volScalarField Q_;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- Calculate Intensity of the laser at p(t, theta) [W/m2]
        scalar calculateIp(scalar r, scalar theta);

        // Normal surface on the interface
        tmp<volVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;


        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;


        //- Disallow default bitwise copy construct
        laserDTRM(const laserDTRM&);

        //- Disallow default bitwise assignment
        void operator=(const laserDTRM&);


public:

    //- Runtime type information
    TypeName("laserDTRM");


    // Constructors

        //- Construct from components
        laserDTRM(const volScalarField& T);

        //- Construct from components
        laserDTRM(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~laserDTRM();


    // Member functions


        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<DimensionedField<scalar, volMesh> > Ru() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
