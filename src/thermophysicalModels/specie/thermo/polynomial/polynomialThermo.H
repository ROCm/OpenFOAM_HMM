/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::polynomialThermo

Description
    Thermodynamics package templated on the equation of state, using polynomial
    functions for cp and h

SourceFiles
    polynomialThermoI.H
    polynomialThermo.C

\*---------------------------------------------------------------------------*/

#ifndef polynomialThermo_H
#define polynomialThermo_H

#include "scalar.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState, int PolySize>
class polynomialThermo;

template<class EquationOfState, int PolySize>
inline polynomialThermo<EquationOfState, PolySize> operator+
(
    const polynomialThermo<EquationOfState, PolySize>&,
    const polynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline polynomialThermo<EquationOfState, PolySize> operator-
(
    const polynomialThermo<EquationOfState, PolySize>&,
    const polynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline polynomialThermo<EquationOfState, PolySize> operator*
(
    const scalar,
    const polynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
inline polynomialThermo<EquationOfState, PolySize> operator==
(
    const polynomialThermo<EquationOfState, PolySize>&,
    const polynomialThermo<EquationOfState, PolySize>&
);

template<class EquationOfState, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const polynomialThermo<EquationOfState, PolySize>&
);


/*---------------------------------------------------------------------------*\
                      Class polynomialThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState, int PolySize>
class polynomialThermo
:
    public EquationOfState
{
    // Private data

        //- Heat of formation [J/kg]
        scalar Hf_;

        //- Standard entropy [J/(kg.K)]
        scalar Sf_;

        //- Specific heat at constant pressure
        Polynomial<PolySize> cpPolynomial_;

        //- Enthalpy - derived from cp
        typename Polynomial<PolySize>::intPolyType dhPolynomial_;

        //- Entropy - derived from cp
        Polynomial<PolySize> sPolynomial_;


    // Private member functions

        //- Construct from components
        inline polynomialThermo
        (
            const EquationOfState& pt,
            const scalar Hf,
            const scalar Sf,
            const Polynomial<PolySize>& cpPoly,
            const typename Polynomial<PolySize>::intPolyType& hPoly,
            const Polynomial<PolySize>& sPoly
        );


public:

    // Constructors

        //- Construct from dictionary
        polynomialThermo(Istream& is);

        //- Construct as a named copy
        inline polynomialThermo(const word&, const polynomialThermo&);


    // Member Functions

        //- Heat capacity at constant pressure [J/(kmol K)]
        inline scalar cp(const scalar T) const;

        //- Enthalpy [J/kmol]
        inline scalar h(const scalar T) const;

        //- Sensible enthalpy [J/kmol]
        inline scalar hs(const scalar T) const;

        //- Chemical enthalpy [J/kmol]
        inline scalar hc() const;

        //- Entropy [J/(kmol K)]
        inline scalar s(const scalar T) const;


    // Member operators

        inline void operator+=(const polynomialThermo&);
        inline void operator-=(const polynomialThermo&);


    // Friend operators

        friend polynomialThermo operator+ <EquationOfState, PolySize>
        (
            const polynomialThermo&,
            const polynomialThermo&
        );

        friend polynomialThermo operator- <EquationOfState, PolySize>
        (
            const polynomialThermo&,
            const polynomialThermo&
        );

        friend polynomialThermo operator* <EquationOfState, PolySize>
        (
            const scalar,
            const polynomialThermo&
        );

        friend polynomialThermo operator== <EquationOfState, PolySize>
        (
            const polynomialThermo&,
            const polynomialThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState, PolySize>
        (
            Ostream&,
            const polynomialThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polynomialThermoI.H"

#ifdef NoRepository
#   include "polynomialThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
