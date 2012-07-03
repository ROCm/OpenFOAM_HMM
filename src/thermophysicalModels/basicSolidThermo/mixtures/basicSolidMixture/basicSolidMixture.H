/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::basicSolidMixture

Description
    Foam::basicSolidMixture

\*---------------------------------------------------------------------------*/

#ifndef basicSolidMixture_H
#define basicSolidMixture_H

#include "volFields.H"
#include "speciesTable.H"
#include "PtrList.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class basicSolidMixture Declaration
\*---------------------------------------------------------------------------*/

class basicSolidMixture
{

protected:

    // Protected data

        typedef speciesTable solidsTable;

        //- The names of the solids
        solidsTable components_;

        //- Solid mass fractions
        PtrList<volScalarField> Y_;


public:

    //- The base class of the mixture
    typedef basicSolidMixture basicMixtureType;


    // Constructors

        //- Construct from word list and mesh
        basicSolidMixture
        (
            const wordList& solidNames,
            const fvMesh&
        );


    //- Destructor
    virtual ~basicSolidMixture()
    {}


    // Member Functions

        //- Return the solid table
        const solidsTable& components() const
        {
            return components_;
        }

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;


        // Per specie thermo properties

            //- Heat capacity at constant pressure [J/(kg K)]
            virtual scalar Cp
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant volume [J/(kg K)]
            virtual scalar Cv
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Chemical enthalpy [J/kg]
            virtual scalar Hc(const label specieI) const = 0;

            //- Density [Kg/m3]
            virtual scalar rho
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const = 0;


       // Per specie transport properties

                //- Thermal conductivity [W/m/K]
                virtual scalar kappa
                (
                    const label specieI,
                    const scalar p,
                    const scalar T
                ) const = 0;

                //- An-isotropic thermal conductivity [W/m/K]
                virtual vector Kappa
                (
                    const label specieI,
                    const scalar p,
                    const scalar T
                ) const = 0;

                //- Thermal diffusivity of enthalpy/internal energy [kg/m/s]
                virtual scalar alpha
                (
                    const label specieI,
                    const scalar p,
                    const scalar T
                ) const = 0;


        // Cell based properties with on the run mixture calculation

            //- Density
            virtual scalar rho(scalar p, scalar T, label celli) const = 0;

            //- Absorption coefficient
            virtual scalar kappaRad(scalar p, scalar T, label celli) const = 0;

            //- Scatter coefficient
            virtual scalar sigmaS(scalar p, scalar T, label celli) const = 0;

            //- Thermal conductivity
            virtual scalar kappa(scalar p, scalar T, label celli) const = 0;

            //- Emissivity coefficient
            virtual scalar emissivity
            (
                scalar p, scalar T, label celli
            ) const = 0;

            //- Specific heat capacity
            virtual scalar Cp(scalar p, scalar T, label celli) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "basicSolidMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
