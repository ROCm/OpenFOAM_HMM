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
    Foam::constSolidThermo

Description
    Constant properties thermodynamics package for solids


SourceFiles
    constSolidThermoI.H
    constSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef constSolidThermo_H
#define constSolidThermo_H

#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class rhoType> class constSolidThermo;

template<class rhoType>
inline constSolidThermo<rhoType> operator*
(
    const scalar,
    const constSolidThermo<rhoType>&
);

/*
template<class rhoType>
inline constSolidThermo<rhoType> operator+
(
    const constSolidThermo<rhoType>&,
    const constSolidThermo<rhoType>&
);

template<class rhoType>
inline constSolidThermo<rhoType> operator-
(
    const constSolidThermo<rhoType>&,
    const constSolidThermo<rhoType>&
);



template<class rhoType>
inline constSolidThermo<rhoType> operator==
(
    const constSolidThermo<rhoType>&,
    const constSolidThermo<rhoType>&
);
*/
template<class rhoType>
Ostream& operator<<
(
    Ostream&,
    const constSolidThermo<rhoType>&
);


/*---------------------------------------------------------------------------*\
                      Class constSolidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class rhoType>
class constSolidThermo
:
    public rhoType
{
    // Private data

        scalar Cp_;
        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline constSolidThermo
        (
            const rhoType& st,
            const scalar Cp,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from dictionary
        constSolidThermo(const dictionary&);

        //- Construct as named copy
        inline constSolidThermo(const word&, const constSolidThermo&);

        //- Construct and return a clone
        inline autoPtr<constSolidThermo> clone() const;

        //- Selector from Istream
        inline static autoPtr<constSolidThermo> New(dictionary& is);


    // Member Functions

         //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;

        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar p, const scalar T) const;

            //- Absolute enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;


    // Member operators

        inline constSolidThermo& operator=(const constSolidThermo&);
        inline void operator+=(const constSolidThermo&);
        inline void operator-=(const constSolidThermo&);

    // Friend operators

        friend constSolidThermo operator* <rhoType>
        (
            const scalar,
            const constSolidThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <rhoType>
        (
            Ostream&,
            const constSolidThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "constSolidThermoI.H"
#   include "constSolidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
