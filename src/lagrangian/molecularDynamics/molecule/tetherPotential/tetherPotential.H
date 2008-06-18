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
    Foam::tetherPotential

Description

SourceFiles
    tetherPotentialI.H
    tetherPotential.C

\*---------------------------------------------------------------------------*/

#ifndef tetherPotential_H
#define tetherPotential_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class tetherPotential Declaration
\*---------------------------------------------------------------------------*/

class tetherPotential
{
    // Private data

        // At the moment this is hard coded to be a "harmonicSpring"
        // potential.  In future use templated classes virtual functions,
        // function pointers and all kinds of good stuff to make this a
        // generic tether force,

        word tetherPotentialName_;

        word tetherPotentialType_;

        scalar springConstant_;


public:

    // Constructors

        //- Construct null
        tetherPotential();

        //- Construct from components
        tetherPotential
        (
            const word& tetherPotentialName,
            const word& tetherPotentialType,
            const scalar springConstant_
        );


    // Destructor

        virtual ~tetherPotential();


    // Member Functions

        // Access

            inline const word& tetherPotentialName() const;

            inline const word& tetherPotentialType() const;

            inline scalar springConstant() const;

            scalar force(const scalar rITMag) const;

            scalar energy(const scalar rITMag) const;


        // Write

            virtual void write(Ostream&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const tetherPotential&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetherPotentialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
