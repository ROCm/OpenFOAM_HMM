/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    Foam::WallInteractionSite

Description
    Stores the patch ID and templated data to represent a collision
    with a wall to be passed to the wall model.

SourceFiles
    WallInteractionSiteI.H
    WallInteractionSite.C
    WallInteractionSiteIO.C

\*---------------------------------------------------------------------------*/

#ifndef WallInteractionSite_H
#define WallInteractionSite_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class WallInteractionSite;

template<class Type>
Istream& operator>>(Istream&, WallInteractionSite<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const WallInteractionSite<Type>&);


/*---------------------------------------------------------------------------*\
                         Class WallInteractionSite Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class WallInteractionSite
{
    // Private data

        //- Index of originating patch
        label patchI_;

        //- Wall data
        Type wallData_;


public:

    // Constructors

        //- Construct null
        WallInteractionSite();

        //- Construct from components
        WallInteractionSite
        (
            label patchI,
            const Type& wallData
        );


    //- Destructor
    ~WallInteractionSite();


    // Member Functions

            //- Return access to the patch index
            inline label patchIndex() const;

            //- Return non-const access to the patch index
            inline label& patchIndex();

            //- Return access to wall data
            inline const Type& wallData() const;

            //- Return non-const access to wall data
            inline Type& wallData();


    // Member Operators

        bool operator==(const WallInteractionSite<Type>&) const;
        bool operator!=(const WallInteractionSite<Type>&) const;


    // IOstream Operators

        friend Istream& operator>> <Type>
        (Istream&, WallInteractionSite<Type>&);

        friend Ostream& operator<< <Type>
        (Ostream&, const WallInteractionSite<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "WallInteractionSiteI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WallInteractionSite.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
