/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Foam::pairPotential

Description

SourceFiles
    pairPotentialI.H
    pairPotential.C

\*---------------------------------------------------------------------------*/

#ifndef pairPotential_H
#define pairPotential_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                        Class pairPotential Declaration
\*---------------------------------------------------------------------------*/

class pairPotential
{
    // Private data

        // At the moment this is hard coded to be a "shiftedLennardJones"
        // potential.  In future use templated classes virtual functions,
        // function pointers and all kinds of good stuff to make this a
        // generic pair force,

        // Not strictly necessary, but useful to keep track of what's what.
        word pairPotentialName_;

        // To be used when runTime selection of potential function is working.
        word pairPotentialType_;

        scalar sigma_;

        scalar epsilon_;

        scalar rCut_;

        scalar rCutSqr_;

        // rCutShiftedForcePotentialConstTerm
        scalar rCutSFPotConst_;

        // rCutShiftedForceConstTerm
        scalar rCutSFConst_;


public:

    // Constructors

        //- Construct null
        pairPotential();

        //- Construct from components
        pairPotential
        (
            const word& pairPotentialName,
            const word& pairPotentialType,
            const scalar sigma,
            const scalar epsilon,
            const scalar rCut
        );


    // Destructor

        virtual ~pairPotential();


    // Member Functions

        // Access

            inline const word& pairPotentialName() const;

            inline const word& pairPotentialType() const;

            inline scalar sigma() const;

            inline scalar epsilon() const;

            inline scalar rCut() const;

            inline scalar rCutSqr() const;


        // Write

            virtual void write(Ostream&) const;

            scalar force(const scalar rIJMag) const;

            scalar energy(const scalar rIJMag) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const pairPotential&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pairPotentialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
