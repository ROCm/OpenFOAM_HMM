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
    Foam::NoInjection

Description
    Place holder for 'none' option

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
protected:

    // Protected member functions

        //- Number of parcels to introduce over the time step
        label parcelsToInject
        (
            const scalar,
            const scalar
        ) const;

        //- Volume of parcels to introduce over the time step
        scalar volumeToInject
        (
            const scalar,
            const scalar
        ) const;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        NoInjection
        (
            const dictionary&,
            CloudType&
        );


    //- Destructor
    virtual ~NoInjection();


    // Member Functions

        //- Flag to indicate whether model activates injection model
        bool active() const;

        //- Return the end-of-injection time
        scalar timeEnd() const;


        // Injection geometry

            //- Set the injection position and owner cell
            void setPositionAndCell
            (
                const label iParcel,
                const scalar time,
                const polyMeshInfo& meshInfo,
                vector& position,
                label& cellOwner
            );

            //- Return the velocity of the parcel to introduce at a time
            vector velocity
            (
                const label,
                const scalar time,
                const polyMeshInfo& meshInfo
            );

            //- Return the diameter of the parcel to introduce at a time
            scalar d0
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
