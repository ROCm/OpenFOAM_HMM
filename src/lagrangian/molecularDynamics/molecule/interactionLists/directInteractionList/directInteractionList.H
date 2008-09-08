/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::directInteractionList

Description

SourceFiles
    directInteractionListI.H
    directInteractionList.C

\*---------------------------------------------------------------------------*/

#ifndef directInteractionList_H
#define directInteractionList_H

#include "interactionLists.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class directInteractionList Declaration
\*---------------------------------------------------------------------------*/

class directInteractionList
:
    public labelListList
{
    // Private data

        const interactionLists& interactionLists_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        directInteractionList(const directInteractionList&);

        //- Disallow default bitwise assignment
        void operator=(const directInteractionList&);

public:

    // Constructors

        //- Construct new data
        directInteractionList
        (
            const interactionLists& interactionLists,
            const scalar& rCutMax
        );

        //- Construct from file
        directInteractionList
        (
            const interactionLists& interactionLists,
            const IOobject& io
        );

    // Destructor

        ~directInteractionList();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // IOstream Operators

        friend Istream& operator>>(Istream&, directInteractionList&);
        friend Ostream& operator<<(Ostream&, const directInteractionList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "directInteractionListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
