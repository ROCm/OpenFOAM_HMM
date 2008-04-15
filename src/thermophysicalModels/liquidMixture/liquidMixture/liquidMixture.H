/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Foam::liquidMixture

Description
    A mixture of liquids.

    The theory in this class is to a very large extent based on the book:
    'The Properties of Gases & Liquids' 4th ed.
    by Reid, Prausnitz and Poling

    For now it does not do much, since the perfect gas equation is used.

\*---------------------------------------------------------------------------*/

#ifndef liquidMixture_H
#define liquidMixture_H

#include "word.H"
#include "scalarField.H"
#include "PtrList.H"
#include "liquid.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class liquidMixture Declaration
\*---------------------------------------------------------------------------*/

class liquidMixture
{
    // Private data

        // maximum reduced temperature
        static const scalar TrMax;

        //- The names of the liquids
        List<word> components_;

        //- The liquid properties
        PtrList<liquid> properties_;


public:


    // Constructors

        //- Construct from dictionary
        liquidMixture
        (
            const dictionary&
        );


    // Selectors

        //- Select construct from dictionary
        static autoPtr<liquidMixture> New
        (
            const dictionary&
        );


    // Member Functions

        inline const List<word>& components() const
        {
            return components_;
        }

        //- Return the liquid properties
        inline const PtrList<liquid>& properties() const
        {
            return properties_;
        }

        //- Calculate the critical temperature of mixture
        scalar Tc(const scalarField& x) const;

        //- Return pseudocritical temperature according to Kay's rule
        scalar Tpc(const scalarField& x) const;

        //- Return pseudocritical pressure (modified Prausnitz and Gunn)
        scalar Ppc(const scalarField& x) const;

        //- Return mixture accentric factor
        scalar omega(const scalarField& x) const;

        //- Return the surface molar fractions
        scalarField Xs
        (
            const scalar p,
            const scalar Tg,
            const scalar Tl,
            const scalarField& xg,
            const scalarField& xl
        ) const;


        //- Calculate the mean molecular weigth [kg/kmol]
        // from mole fractions
        scalar W(const scalarField& x) const;

        //- Returns the mass fractions, given mole fractions
        scalarField Y(const scalarField& X) const;

        //- Returns the mole fractions, given mass fractions
        scalarField X(const scalarField& Y) const;

        //- Calculate the mixture density [kg/m^3]
        scalar rho
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture vapour pressure [Pa]
        scalar pv
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture latent heat [J/kg]
        scalar hl
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture heat capacity [J/(kg K)]
        scalar cp
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Estimate mixture surface tension [N/m]
        scalar sigma
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Calculate the mixture viscosity [Pa s]
        scalar mu
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        //- Estimate thermal conductivity  [W/(m K)]
        // Li's method, Eq. 10-12.27 - 10.12-19
        scalar K
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;

        scalar D
        (
            const scalar p,
            const scalar T,
            const scalarField& x
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
