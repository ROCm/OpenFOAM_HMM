/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::multiComponentSolidMixture

Description
    Foam::multiComponentSolidMixture

\*---------------------------------------------------------------------------*/

#ifndef multiComponentSolidMixture_H
#define multiComponentSolidMixture_H

#include "PtrList.H"
#include "autoPtr.H"
#include "basicSolidMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class multiComponentSolidMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class multiComponentSolidMixture
:
    public basicSolidMixture
{

    // Private data

        //- Solid data
        PtrList<ThermoType> solidData_;

        //- Temporary storage for the mass cell/face mixture thermo data
        mutable ThermoType mixture_;

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        mutable ThermoType mixtureVol_;


    // Private Member Functions

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();

        //- Return molar fraction for specieI in celli
        scalar X(label specieI, label celli, scalar p, scalar T) const;

        //- Construct solid data
        const ThermoType& constructSolidData(const dictionary&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary and mesh
        multiComponentSolidMixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~multiComponentSolidMixture()
    {}


    // Member Functions


        const ThermoType& cellMixture(const label celli) const;

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        //- Return the raw solid data
        const PtrList<ThermoType>& solidData() const
        {
            return solidData_;
        }

        //- Read dictionary
        void read(const dictionary&);


         // Per specie thermo properties

            //- Heat capacity at constant pressure [J/(kg K)]
            virtual scalar Cp
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const;

            //- Heat capacity at constant volume [J/(kg K)]
            virtual scalar Cv
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const;

            //- Chemical enthalpy [J/kg]
            virtual scalar Hc(const label specieI) const;

            //- Density [Kg/m3]
            virtual scalar rho
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const;


            // Per specie transport properties

                //- Thermal conductivity [W/m/K]
                virtual scalar kappa
                (
                    const label specieI,
                    const scalar p,
                    const scalar T
                ) const;

                //- An-isotropic thermal conductivity [W/m/K]
                virtual vector Kappa
                (
                    const label specieI,
                    const scalar p,
                    const scalar T
                ) const;

                //- Thermal diffusivity of enthalpy/internal energy [kg/m/s]
                virtual scalar alpha
                (
                    const label specieI,
                    const scalar p,
                    const scalar T
                ) const;


        // Cell based properties with on the run mixture calculation.

            //- Density
            virtual scalar rho(scalar p, scalar T, label celli) const;

            //- Thermal conductivity
            virtual scalar kappa(scalar p, scalar T, label celli) const;

            //- Specific heat capacity
            virtual scalar Cp(scalar p, scalar T, label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "multiComponentSolidMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
