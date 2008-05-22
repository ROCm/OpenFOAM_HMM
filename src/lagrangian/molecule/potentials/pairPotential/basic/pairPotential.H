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

    At the moment this is hard coded to be a shifted force
    "maitlandSmith" potential.
    In future use templated classes virtual functions,
    function pointers and all kinds of good stuff to make
    this a generic pair force,

SourceFiles
    pairPotentialI.H
    pairPotential.C

\*---------------------------------------------------------------------------*/

#ifndef pairPotential_H
#define pairPotential_H

#include "vector.H"
#include "dictionary.H"
#include "List.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pairPotential Declaration
\*---------------------------------------------------------------------------*/

class pairPotential
{
    // Private data

        List<scalar> forceLookup_;

        List<scalar> energyLookup_;

        scalar m_;

        scalar gamma_;

        scalar rm_;

        scalar epsilon_;

        scalar rCut_;

        scalar rCutSqr_;

        scalar u_at_rCut_;

        scalar du_by_dr_at_rCut_;

        scalar rMin_;

        scalar dr_;


    // Private Member Functions

        void setLookupTables();

        void setConstants();


public:

    // Constructors

        //- Construct null
        pairPotential();

        //- Construct from dictionary
        pairPotential(const dictionary& pPDict);

        //- Construct from components
        pairPotential
        (
            const scalar m,
            const scalar gamma,
            const scalar rm,
            const scalar epsilon,
            const scalar rCut,
            const scalar rMin,
            const scalar dr
        );


    // Destructor

        virtual ~pairPotential();


    // Member Functions

        // Access

            inline const scalar m() const;

            inline const scalar gamma() const ;

            inline const scalar rm() const;

            inline const scalar epsilon() const;

            inline const scalar rCut() const;

            inline const scalar rCutSqr() const;

            inline const scalar rMin() const;

            inline const scalar dr() const;


        // Write

            virtual void write(Ostream&) const;

            scalar n(const scalar rIJMag) const;

            scalar force(const scalar rIJMag) const;

            scalar forceLookup(const scalar rIJMag) const;

            List<Pair<scalar> > forceTable() const;

            scalar energy(const scalar rIJMag) const;

            scalar energyLookup(const scalar rIJMag) const;

            List<Pair<scalar> > energyTable() const;


    // Friend Operators

        inline friend bool operator==
        (
            const pairPotential& a,
            const pairPotential& b
        );

        inline friend bool operator!=
        (
            const pairPotential& a,
            const pairPotential& b
        );


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
