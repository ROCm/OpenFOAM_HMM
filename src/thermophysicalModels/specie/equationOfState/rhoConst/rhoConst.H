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
    Foam::rhoConst

Description
    RhoConst (rho = const) of state.

SourceFiles
    rhoConstI.H
    rhoConst.C

\*---------------------------------------------------------------------------*/

#ifndef rhoConst_H
#define rhoConst_H

#include "specie.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class rhoConst Declaration
\*---------------------------------------------------------------------------*/

class rhoConst
:
    public specie
{
    // Private data

        //- Density
        scalar rho_;


public:

    // Constructors

        //- Construct from components
        inline rhoConst(const specie& sp, const scalar rho);

        //- Construct from Istream
        rhoConst(Istream&);

        //- Construct from dictionary
        rhoConst(const dictionary& dict);

        //- Construct as named copy
        inline rhoConst(const word& name, const rhoConst&);

        //- Construct and return a clone
        inline autoPtr<rhoConst> clone() const;

        // Selector from Istream
        inline static autoPtr<rhoConst> New(Istream& is);


    // Member functions

        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = true;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const rhoConst&);
        inline void operator-=(const rhoConst&);

        inline void operator*=(const scalar);


    // Friend operators

        inline friend rhoConst operator+
        (
            const rhoConst&,
            const rhoConst&
        );

        inline friend rhoConst operator-
        (
            const rhoConst&,
            const rhoConst&
        );

        inline friend rhoConst operator*
        (
            const scalar s,
            const rhoConst&
        );

        inline friend rhoConst operator==
        (
            const rhoConst&,
            const rhoConst&
        );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const rhoConst&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rhoConstI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
