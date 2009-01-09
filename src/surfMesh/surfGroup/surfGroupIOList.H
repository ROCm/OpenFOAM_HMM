/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::surfGroupIOList

Description
    IOobject for a surfGroupList

SourceFiles
    surfGroupIOList.C

\*---------------------------------------------------------------------------*/

#ifndef surfGroupIOList_H
#define surfGroupIOList_H

#include "surfGroupList.H"
#include "regIOobject.H"
#include "faceList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class surfGroupIOList Declaration
\*---------------------------------------------------------------------------*/

class surfGroupIOList
:
    public surfGroupList,
    public regIOobject
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        surfGroupIOList(const surfGroupIOList&);

        //- Disallow default bitwise assignment
        void operator=(const surfGroupIOList&);


public:

    //- Runtime type information
    TypeName("surfGroupIOList");


    // Constructors

        //- Construct from IOobject
        explicit surfGroupIOList(const IOobject& io);

        //- Construct from IOobject
        surfGroupIOList(const IOobject& io, const surfGroupList&);

    // Destructor

        ~surfGroupIOList();


    // Member Functions

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const surfGroupIOList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
