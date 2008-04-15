/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Foam::NoInjection

Description
    Place holder for 'no injection' option

SourceFiles
    NoInjection.C

\*---------------------------------------------------------------------------*/

#ifndef noInjection_H
#define noInjection_H

#include "InjectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class NoInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoInjection
:
    public InjectionModel<CloudType>
{

public:

    //- Runtime type information
    TypeName("NoInjection");


    // Constructors

        //- Construct from components
        NoInjection
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Destructor

        ~NoInjection();


    // Member Functions

        bool active() const;

        scalar timeStart() const;

        scalar timeEnd() const;

        vector position
        (
            const label,
            const scalar,
            const polyMeshInfo&,
            Random&
        ) const;

        label nParcelsToInject
        (
            const label,
            const scalar,
            const scalar
        ) const;

        scalar volume
        (
            const scalar,
            const scalar,
            const polyMeshInfo&
        ) const;

        scalar volumeFraction
        (
            const scalar,
            const scalar
        ) const;

        scalar d0
        (
            const label,
            const scalar
        ) const;

        vector velocity
        (
            const label,
            const scalar
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
