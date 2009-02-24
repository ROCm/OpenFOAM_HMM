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
    Foam::DsmcCloud

Description
    Templated base class for dsmc cloud

SourceFiles
    DsmcCloudI.H
    DsmcCloud.C

\*---------------------------------------------------------------------------*/

#ifndef DsmcCloud_H
#define DsmcCloud_H

#include "Cloud.H"
#include "DsmcBaseCloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"

#include "IntegrationSchemesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class CollisionModel;

template<class CloudType>
class InjectionModel;

template<class CloudType>
class WallInteractionModel;

/*---------------------------------------------------------------------------*\
                       Class DsmcCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class DsmcCloud
:
    public Cloud<ParcelType>,
    public DsmcBaseCloud
{

private:

    // Private data

        //- Cloud type - used to set the name of the parcel properties
        //  dictionary by appending "Properties"
        const word cloudType_;

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- A data structure holding which particles are in which cell
        List<DynamicList<ParcelType*> > cellOccupancy_;

        //- Parcel constant properties - one for each species
        List<typename ParcelType::constantProperties> constProps_;

        //- Random number generator
        Random rndGen_;

        //- Interpolation schemes dictionary
        dictionary interpolationSchemes_;


        // References to the cloud sub-models

            //- Injector model
            autoPtr<CollisionModel<DsmcCloud<ParcelType> > >
                collisionModel_;

            //- Injector model
            autoPtr<InjectionModel<DsmcCloud<ParcelType> > >
                injectionModel_;

            //- Wall interaction model
            autoPtr<WallInteractionModel<DsmcCloud<ParcelType> > >
                wallInteractionModel_;


        // Reference to the particle integration schemes

            //- Velocity integration
            autoPtr<vectorIntegrationScheme> UIntegrator_;


    // Private Member Functions

        //- Record which particles are in which cell
        void buildCellOccupancy()

        //- Disallow default bitwise copy construct
        DsmcCloud(const DsmcCloud&);

        //- Disallow default bitwise assignment
        void operator=(const DsmcCloud&);


public:

    // Constructors

        //- Construct given carrier gas fields
        DsmcCloud
        (
            const word& cloudType,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~DsmcCloud();


    // Member Functions

        // Access

            // References to the mesh and databases

                //- Return the cloud type
                inline const word& cloudType() const;

                //- Return refernce to the mesh
                inline const fvMesh& mesh() const;

                //- Return particle properties dictionary
                inline const IOdictionary& particleProperties() const;

            //- Return refernce to the random object
            inline Random& rndGen();


            // Sub-models

                //- Return reference to injection model
                inline const CollisionModel<DsmcCloud<ParcelType> >&
                    collision() const;

                //- Return reference to injection model
                inline const InjectionModel<DsmcCloud<ParcelType> >&
                    injection() const;

                inline InjectionModel<DsmcCloud<ParcelType> >&
                    injection();

                //- Return reference to wall interaction model
                inline const WallInteractionModel<DsmcCloud<ParcelType> >&
                    wallInteraction() const;


            // Integration schemes

                //-Return reference to velocity integration
                inline const vectorIntegrationScheme& UIntegrator() const;



        // Check

            //- Total mass injected
            inline scalar massInjected() const;

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Print cloud information
            void info() const;

            //- Dump particle positions to .obj file
            void dumpParticlePositions() const;


            // Fields

                //- Return the real particle number density field
                inline const tmp<volScalarField> rhoN() const;

                //- Return the particle mass density field
                inline const tmp<volScalarField> rhoM() const;

                //- Return the field of number of DSMC particles
                inline const tmp<volScalarField> rhoNP() const;

                //- Return the velocity field
                inline const tmp<volVectorField> U() const;

                //- Return the temperature field
                inline const tmp<volScalarField> T() const;


        // Cloud evolution functions

            //- Add new parcel
            void addNewParcel
            (
                const vector& position,
                const vector& U,
                const label cellId,
                const label speciesId,
            );

            //- Evolve the cloud (move, collide)
            void evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DsmcCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DsmcCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
