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
    Foam::interactionLists

Description

SourceFiles
    interactionListsI.H
    interactionLists.C
    interactionListsIO.C

\*---------------------------------------------------------------------------*/

#ifndef interactionLists_H
#define interactionLists_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class interactionLists;
Istream& operator>>(Istream&, interactionLists&);
Ostream& operator<<(Ostream&, const interactionLists&);


/*---------------------------------------------------------------------------*\
                           Class interactionLists Declaration
\*---------------------------------------------------------------------------*/

class interactionLists
:
    public baseClassName
{
    // Private data

        dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        interactionLists(const interactionLists&);

        //- Disallow default bitwise assignment
        void operator=(const interactionLists&);


public:

    // Static data members

        //- Static data someStaticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        interactionLists();

        //- Construct from components
        interactionLists(const dataType& data);

        //- Construct from Istream
        interactionLists(Istream&);

        //- Construct as copy
        interactionLists(const interactionLists&);


    // Selectors

        //- Select null constructed
        static autoPtr<interactionLists> New();


    // Destructor

        ~interactionLists();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const interactionLists&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, interactionLists&);
        friend Ostream& operator<<(Ostream&, const interactionLists&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interactionListsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
