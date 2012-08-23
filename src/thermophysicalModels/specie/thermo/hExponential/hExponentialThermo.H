/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::hExponentialThermo

Description
    Exponential properties thermodynamics package
    templated into the equationOfState.

SourceFiles
    hExponentialThermoI.H
    hExponentialThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hExponentialThermo_H
#define hExponentialThermo_H

#include "scalar.H"
#include "dictionary.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class equationOfState> class hExponentialThermo;

template<class equationOfState>
inline hExponentialThermo<equationOfState> operator*
(
    const scalar,
    const hExponentialThermo<equationOfState>&
);


template<class equationOfState>
Ostream& operator<<
(
    Ostream&,
    const hExponentialThermo<equationOfState>&
);


/*---------------------------------------------------------------------------*\
                         Class hExponentialThermo Declaration
\*---------------------------------------------------------------------------*/

template<class equationOfState>
class hExponentialThermo
:
    public equationOfState
{

    // Private data

        scalar c0_;
        scalar n0_;
        scalar Tref_;
        scalar Hf_;


    // Private Member Functions

        //- Check given temperature is within the range of the fitted coeffs
        inline void checkT(const scalar T) const;

        //- Integrate Cp expression
        inline scalar integrateCp(const scalar T) const;


public:

    // Constructors

        //- Construct from components
        inline hExponentialThermo
        (
            const equationOfState& st,
            const scalar c0,
            const scalar n0,
            const scalar Tref,
            const scalar Hf
        );

        //- Construct from dictionary
        hExponentialThermo(const dictionary&);

        //- Construct as a named copy
        inline hExponentialThermo
        (
            const word&,
            const hExponentialThermo&
        );


    // Member Functions

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar cp(const scalar p, const scalar T) const;

            //- Absolute enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar hs(const scalar p, const scalar T) const;

            //- chemical enthalpy [J/kg]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;


    // Member operators

        inline hExponentialThermo& operator=
        (
            const hExponentialThermo&
        );
        inline void operator+=(const hExponentialThermo&);
        inline void operator-=(const hExponentialThermo&);


    // Friend operators


        friend hExponentialThermo operator* <equationOfState>
        (
            const scalar,
            const hExponentialThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <equationOfState>
        (
            Ostream&,
            const hExponentialThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "hExponentialThermoI.H"
#   include "hExponentialThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
