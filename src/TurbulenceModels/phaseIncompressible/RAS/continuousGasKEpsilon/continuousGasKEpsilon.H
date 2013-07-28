/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::RASModels::continuousGasKEpsilon

Group
    grpRASTurbulence

Description

SourceFiles
    continuousGasKEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef continuousGasKEpsilon_H
#define continuousGasKEpsilon_H

#include "kEpsilon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class continuousGasKEpsilon Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class continuousGasKEpsilon
:
    public kEpsilon<BasicTurbulenceModel>
{
    // Private data

        mutable const turbulenceModel *liquidTurbulencePtr_;

        volScalarField nutEff_;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar alphaInversion_;


    // Protected member functions

        virtual void correctNut();
        tmp<volScalarField> phaseTransferCoeff() const;
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("continuousGasKEpsilon");


    // Constructors

        //- Construct from components
        continuousGasKEpsilon
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~continuousGasKEpsilon()
    {}


    // Member Functions

        //- Return the turbulence model for the liquid phase
        const turbulenceModel& liquidTurbulence() const;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const;

        //- Return the effective density for the stress
        virtual tmp<volScalarField> rhoEff() const;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "continuousGasKEpsilon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
