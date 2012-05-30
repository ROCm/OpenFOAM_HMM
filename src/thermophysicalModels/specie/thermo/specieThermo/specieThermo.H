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
    Foam::specieThermo

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type thermo.  All other
    properties are derived from these primitive functions.

SourceFiles
    specieThermoI.H
    specieThermo.C

\*---------------------------------------------------------------------------*/

#ifndef specieThermo_H
#define specieThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo, template<class> class Type> class specieThermo;

template<class Thermo, template<class> class Type>
inline specieThermo<Thermo, Type> operator+
(
    const specieThermo<Thermo, Type>&,
    const specieThermo<Thermo, Type>&
);

template<class Thermo, template<class> class Type>
inline specieThermo<Thermo, Type> operator-
(
    const specieThermo<Thermo, Type>&,
    const specieThermo<Thermo, Type>&
);

template<class Thermo, template<class> class Type>
inline specieThermo<Thermo, Type> operator*
(
    const scalar,
    const specieThermo<Thermo, Type>&
);

template<class Thermo, template<class> class Type>
inline specieThermo<Thermo, Type> operator==
(
    const specieThermo<Thermo, Type>&,
    const specieThermo<Thermo, Type>&
);

template<class Thermo, template<class> class Type>
Ostream& operator<<
(
    Ostream&,
    const specieThermo<Thermo, Type>&
);


/*---------------------------------------------------------------------------*\
                           Class specieThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, template<class> class Type>
class specieThermo
:
    public Thermo,
    public Type<specieThermo<Thermo, Type> >
{
    // Private data

        //- Convergence tolerance of energy -> temperature inversion functions
        static const scalar tol_;

        //- Max number of iterations in energy->temperature inversion functions
        static const int maxIter_;


    // Private Member Functions

        //- Return the temperature corresponding to the value of the
        //  thermodynamic property f, given the function f = F(T) and dF(T)/dT
        inline scalar T
        (
            scalar f,
            scalar T0,
            scalar (specieThermo::*F)(const scalar) const,
            scalar (specieThermo::*dFdT)(const scalar) const,
            scalar (specieThermo::*limit)(const scalar) const
        ) const;


public:

    // Constructors

        //- construct from components
        inline specieThermo(const Thermo& sp);

        //- Construct from Istream
        specieThermo(Istream&);

        //- Construct from dictionary
        specieThermo(const dictionary& dict);

        //- Construct as named copy
        inline specieThermo(const word& name, const specieThermo&);


    // Member Functions

        // Fundamental properties
        // (These functions must be provided in derived types)

            // Heat capacity at constant pressure [J/(kmol K)]
            // scalar cp(const scalar) const;

            // Absolute Enthalpy [J/kmol]
            // scalar ha(const scalar) const;

            // Sensible enthalpy [J/kmol]
            // scalar hs(const scalar) const;

            // Chemical enthalpy [J/kmol]
            // scalar hc(const scalar) const;

            // Entropy [J/(kmol K)]
            // scalar s(const scalar) const;


        // Calculate and return derived properties
        // (These functions need not provided in derived types)

            // Mole specific properties

                //- Name of Enthalpy/Internal energy
                static inline word heName();

                //- Enthalpy/Internal energy [J/kmol]
                inline scalar he(const scalar T) const;

                //- Heat capacity at constant volume [J/(kmol K)]
                inline scalar cv(const scalar T) const;

                //- Heat capacity at constant pressure/volume [J/(kmol K)]
                inline scalar cpv(const scalar T) const;

                //- gamma = cp/cv []
                inline scalar gamma(const scalar T) const;

                //- Ratio of heat capacity at constant pressure to that at
                //  constant pressure/volume []
                inline scalar cpBycpv(const scalar T) const;

                //- Sensible internal energy [J/kmol]
                inline scalar es(const scalar T) const;

                //- Absolute internal energy [J/kmol]
                inline scalar ea(const scalar T) const;

                //- Gibbs free energy [J/kmol]
                inline scalar g(const scalar T) const;

                //- Helmholtz free energy [J/kmol]
                inline scalar a(const scalar T) const;


            // Mass specific properties

                //- Heat capacity at constant pressure [J/(kg K)]
                inline scalar Cp(const scalar T) const;

                //- Heat capacity at constant volume [J/(kg K)]
                inline scalar Cv(const scalar T) const;

                //- Heat capacity at constant pressure/volume [J/(kg K)]
                inline scalar Cpv(const scalar T) const;

                //- Enthalpy/Internal energy [J/kg]
                inline scalar HE(const scalar T) const;

                //- Enthalpy [J/kg]
                inline scalar H(const scalar T) const;

                //- Sensible enthalpy [J/kg]
                inline scalar Hs(const scalar T) const;

                //- Chemical enthalpy [J/kg]
                inline scalar Hc() const;

                //- Absolute Enthalpy [J/kg]
                inline scalar Ha(const scalar T) const;

                //- Entropy [J/(kg K)]
                inline scalar S(const scalar T) const;

                //- Internal energy [J/kg]
                inline scalar E(const scalar T) const;

                //- Sensible internal energy [J/kg]
                inline scalar Es(const scalar T) const;

                //- Absolute internal energy [J/kg]
                inline scalar Ea(const scalar T) const;

                //- Gibbs free energy [J/kg]
                inline scalar G(const scalar T) const;

                //- Helmholtz free energy [J/kg]
                inline scalar A(const scalar T) const;


        // Equilibrium reaction thermodynamics

            //- Equilibrium constant [] i.t.o fugacities
            //  = PIi(fi/Pstd)^nui
            inline scalar K(const scalar T) const;

            //- Equilibrium constant [] i.t.o. partial pressures
            //  = PIi(pi/Pstd)^nui
            //  For low pressures (where the gas mixture is near perfect) Kp = K
            inline scalar Kp(const scalar T) const;

            //- Equilibrium constant i.t.o. molar concentration
            //  = PIi(ci/cstd)^nui
            //  For low pressures (where the gas mixture is near perfect)
            //  Kc = Kp(pstd/(RR*T))^nu
            inline scalar Kc(const scalar T) const;

            //- Equilibrium constant [] i.t.o. mole-fractions
            //  For low pressures (where the gas mixture is near perfect)
            //  Kx = Kp(pstd/p)^nui
            inline scalar Kx(const scalar T, const scalar p) const;

            //- Equilibrium constant [] i.t.o. number of moles
            //  For low pressures (where the gas mixture is near perfect)
            //  Kn = Kp(n*pstd/p)^nui where n = number of moles in mixture
            inline scalar Kn
            (
                const scalar T,
                const scalar p,
                const scalar n
            ) const;


        // Energy->temperature  inversion functions

            //- Temperature from enthalpy or internal energy
            //  given an initial temperature T0
            inline scalar THE(const scalar H, const scalar T0) const;

            //- Temperature from sensible enthalpy given an initial T0
            inline scalar THs(const scalar Hs, const scalar T0) const;

            //- Temperature from absolute enthalpy
            //  given an initial temperature T0
            inline scalar THa(const scalar H, const scalar T0) const;

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            inline scalar TEs(const scalar E, const scalar T0) const;

            //- Temperature from absolute internal energy
            //  given an initial temperature T0
            inline scalar TEa(const scalar E, const scalar T0) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const specieThermo&);
        inline void operator-=(const specieThermo&);

        inline void operator*=(const scalar);


    // Friend operators

        friend specieThermo operator+ <Thermo, Type>
        (
            const specieThermo&,
            const specieThermo&
        );

        friend specieThermo operator- <Thermo, Type>
        (
            const specieThermo&,
            const specieThermo&
        );

        friend specieThermo operator* <Thermo, Type>
        (
            const scalar s,
            const specieThermo&
        );

        friend specieThermo operator== <Thermo, Type>
        (
            const specieThermo&,
            const specieThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo, Type>
        (
            Ostream&,
            const specieThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "specieThermoI.H"

#ifdef NoRepository
#   include "specieThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
