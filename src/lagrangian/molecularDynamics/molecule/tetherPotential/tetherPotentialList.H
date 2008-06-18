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
    Foam::tetherPotentialList

Description

SourceFiles
    tetherPotentialListI.H
    tetherPotentialList.C

\*---------------------------------------------------------------------------*/

#ifndef tetherPotentialList_H
#define tetherPotentialList_H

#include "ListOps.H"
#include "tetherPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class tetherPotentialList Declaration
\*---------------------------------------------------------------------------*/

class tetherPotentialList
:
    public List<tetherPotential>
{
    // Private data

        List<label> tetherIds_;


    // Private Member Functions

        label tetherPotentialIndex(const label a);

        //- Disallow default bitwise copy construct
        tetherPotentialList(const tetherPotentialList&);


public:

    // Constructors

        //- Construct null
        tetherPotentialList();

        //- Construct from list of tetherIDs
        tetherPotentialList(const List<label> tetherIds);


    // Destructor

        ~tetherPotentialList();


    // Member Functions

        // Access

            inline const List<label>& tetherIds() const;

            inline label nTetherPotentials() const;

            void setSizeAndTetherIds (const List<label>&);

            void addPotential
            (
                const label a,
                const tetherPotential& pot
            );

            tetherPotential& tetherPotentialFunction (const label a);

            scalar force (const label a, const scalar rITMag);

            scalar energy (const label a, const scalar rITMag);


    // Member Operators

        void operator=(const tetherPotentialList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetherPotentialListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
