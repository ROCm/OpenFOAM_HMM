/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::regionModels::areaSurfaceFilmModels::
        perturbedTemperatureDependentContactAngleForce

Description
    Temperature dependent contact angle force with a stochastic perturbation.

    The contact angle in degrees is specified as a \c Foam::Function1 type,
    to enable the use of, e.g. \c constant, \c polynomial, \c table values
    and the stochastic perturbation obtained from a
    \c Foam::distributionModels::distributionModel.

See also
  - Foam::regionModels::areaSurfaceFilmModels::contactAngleForce
  - areaSurfaceFilmModels::temperatureDependentContactAngleForce
  - Foam::regionModels::areaSurfaceFilmModels::distributionContactAngleForce
  - Foam::Function1Types
  - Foam::distributionModel

SourceFiles
    perturbedTemperatureDependentContactAngleForce.C

\*---------------------------------------------------------------------------*/

#ifndef perturbedTemperatureDependentContactAngleForce_H
#define perturbedTemperatureDependentContactAngleForce_H

#include "contactAngleForce.H"
#include "Function1.H"
#include "distributionModel.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
       Class perturbedTemperatureDependentContactAngleForce Declaration
\*---------------------------------------------------------------------------*/

class perturbedTemperatureDependentContactAngleForce
:
    public contactAngleForce
{
    // Private Data

        //- Contact angle function
        autoPtr<Function1<scalar>> thetaPtr_;

        //- Random number generator
        Random rndGen_;

        //- Parcel size PDF model
        const autoPtr<distributionModel> distribution_;


    // Private Member Functions

        //- No copy construct
        perturbedTemperatureDependentContactAngleForce
        (
            const perturbedTemperatureDependentContactAngleForce&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const perturbedTemperatureDependentContactAngleForce&
        ) = delete;


protected:

        //- Return the contact angle field
        virtual tmp<areaScalarField> theta() const;


public:

    //- Runtime type information
    TypeName("perturbedTemperatureDependentContactAngle");


    // Constructors

        //- Construct from surface film model
        perturbedTemperatureDependentContactAngleForce
        (
            liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~perturbedTemperatureDependentContactAngleForce();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
