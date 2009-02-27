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
    Foam::ReactingMultiphaseCloud

Description
    Templated base class for multiphase reacting cloud
    - Adds to reacting cloud
      - Devolatilisatsion
      - multiphase composition

SourceFiles
    ReactingMultiphaseCloudI.H
    ReactingMultiphaseCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingMultiphaseCloud_H
#define ReactingMultiphaseCloud_H

#include "ReactingCloud.H"
#include "reactingMultiphaseCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class DevolatilisationModel;

/*---------------------------------------------------------------------------*\
                   Class ReactingMultiphaseCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingMultiphaseCloud
:
    public ReactingCloud<ParcelType>,
    public reactingMultiphaseCloud
{
    // Private data

        //- Parcel constant properties
        typename ParcelType::constantProperties constProps_;


        // References to the cloud sub-models

            //- Devolatilisation model
            autoPtr
            <
                DevolatilisationModel<ReactingMultiphaseCloud<ParcelType> >
            >
            devolatilisationModel_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ReactingMultiphaseCloud(const ReactingMultiphaseCloud&);

        //- Disallow default bitwise assignment
        void operator=(const ReactingMultiphaseCloud&);


public:

    // Constructors

        //- Construct given carrier gas fields
        ReactingMultiphaseCloud
        (
            const word& cloudType,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            hCombustionThermo& thermo,
            PtrList<specieReactingProperties>& gases
        );


    //- Destructor
    virtual ~ReactingMultiphaseCloud();


    //- Type of parcel the cloud was instantiated for
    typedef ParcelType parcelType;


    // Member Functions

        // Access

            // Sub-models

                //- Return reference to devolatilisation model
                inline const DevolatilisationModel
                <
                    ReactingMultiphaseCloud<ParcelType>
                >&
                devolatilisation() const;


        // Cloud evolution functions

            //- Add new parcel
            void addNewParcel
            (
                const vector& position,
                const label cellId,
                const scalar d,
                const vector& U,
                const scalar nParticles,
                const scalar lagrangianDt
            );

            //- Reset the spray source terms
            void resetSourceTerms();

            //- Evolve the spray (inject, move)
            void evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingMultiphaseCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReactingMultiphaseCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
