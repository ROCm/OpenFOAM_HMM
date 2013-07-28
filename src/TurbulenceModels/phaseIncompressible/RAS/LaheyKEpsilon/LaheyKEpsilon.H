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
    Foam::RASModels::LaheyKEpsilon

Group
    grpRASTurbulence

Description

SourceFiles
    LaheyKEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef LaheyKEpsilon_H
#define LaheyKEpsilon_H

#include "kEpsilon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class LaheyKEpsilon Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class LaheyKEpsilon
:
    public kEpsilon<BasicTurbulenceModel>
{
    // Private data

        mutable const PhaseIncompressibleTurbulenceModel
        <
            typename BasicTurbulenceModel::transportModel
        > *gasTurbulencePtr_;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar alphaInversion_;
            dimensionedScalar Cp_;
            dimensionedScalar Cmub_;


    // Protected member functions

        virtual void correctNut();
        tmp<volScalarField> bubbleG() const;
        tmp<volScalarField> phaseTransferCoeff() const;
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("LaheyKEpsilon");


    // Constructors

        //- Construct from components
        LaheyKEpsilon
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
    virtual ~LaheyKEpsilon()
    {}


    // Member Functions

        //- Return the turbulence model for the gas phase
        const PhaseIncompressibleTurbulenceModel<transportModel>&
        gasTurbulence() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LaheyKEpsilon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
