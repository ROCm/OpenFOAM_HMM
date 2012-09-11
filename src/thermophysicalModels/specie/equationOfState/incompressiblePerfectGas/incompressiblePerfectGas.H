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
    Foam::incompressiblePerfectGas

Description
    Incompressible gas equation of state using a constant reference pressure in
    the perfect gas equation of state rather than the local pressure so that the
    density only varies with temperature and composition.

SourceFiles
    incompressiblePerfectGasI.H
    incompressiblePerfectGas.C

\*---------------------------------------------------------------------------*/

#ifndef incompressiblePerfectGas_H
#define incompressiblePerfectGas_H

#include "specie.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class incompressiblePerfectGas Declaration
\*---------------------------------------------------------------------------*/

class incompressiblePerfectGas
:
    public specie
{
    // Private data

        //- Reference pressure
        scalar pRef_;


public:

    // Constructors

        //- Construct from components
        inline incompressiblePerfectGas(const specie& sp);

        //- Construct from Istream
        incompressiblePerfectGas(Istream&);

        //- Construct from dictionary
        incompressiblePerfectGas(const dictionary& dict);

        //- Construct as named copy
        inline incompressiblePerfectGas
        (
            const word& name,
            const incompressiblePerfectGas&
        );

        //- Construct and return a clone
        inline autoPtr<incompressiblePerfectGas> clone() const;

        // Selector from Istream
        inline static autoPtr<incompressiblePerfectGas> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<incompressiblePerfectGas> New
        (
            const dictionary& dict
        );


    // Member functions

        // Fundamental properties

            //- Return true if the equation of state is incompressible
            //  i.e. rho != f(p)
            static const bool incompressible = true;

            //- Return true if the equation of state is isochoric
            //  i.e. rho = const
            static const bool isochoric = false;

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

        inline void operator+=(const incompressiblePerfectGas&);
        inline void operator-=(const incompressiblePerfectGas&);

        inline void operator*=(const scalar);


    // Friend operators

        inline friend incompressiblePerfectGas operator+
        (
            const incompressiblePerfectGas&,
            const incompressiblePerfectGas&
        );

        inline friend incompressiblePerfectGas operator-
        (
            const incompressiblePerfectGas&,
            const incompressiblePerfectGas&
        );

        inline friend incompressiblePerfectGas operator*
        (
            const scalar s,
            const incompressiblePerfectGas&
        );

        inline friend incompressiblePerfectGas operator==
        (
            const incompressiblePerfectGas&,
            const incompressiblePerfectGas&
        );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const incompressiblePerfectGas&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "incompressiblePerfectGasI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
