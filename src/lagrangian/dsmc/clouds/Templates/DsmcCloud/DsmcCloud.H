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
#include "scalarIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class BinaryCollisionModel;

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
    // Private data

        //- Cloud type - used to set the name of the parcel properties
        //  dictionary by appending "Properties"
        const word cloudType_;

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- A list of unique instances of molecule types in the simulation.
        // The position of an entry in the list maps to the label identifying
        // the typeId, i.e. where typeIdList_ = (N2 O2 CO2)
        // N2 has typeId label = 0, O2 = 1, CO2 = 2.
        List<word> typeIdList_;

        //- Number of real atoms/molecules represented by a parcel
        scalar nParticle_;

        //- A data structure holding which particles are in which cell
        List<DynamicList<ParcelType*> > cellOccupancy_;

        //- An IOField holding the value of (sigmaT * cR)max for each cell (see
        // Bird p220). Initialised with the parcelsm updated as required, and
        // read in on start/restart.
        volScalarField sigmaTcRMax_;

        //- A field holding the remainder from the previous collision selections
        scalarField collisionSelectionRemainder_;

        //- Parcel constant properties - one for each type
        List<typename ParcelType::constantProperties> constProps_;

        //- Random number generator
        Random rndGen_;

        // References to the cloud sub-models

            //- Binary  Collision model
            autoPtr<BinaryCollisionModel<DsmcCloud<ParcelType> > >
                binaryCollisionModel_;

            //- Wall interaction model
            autoPtr<WallInteractionModel<DsmcCloud<ParcelType> > >
                wallInteractionModel_;


    // Private Member Functions

        //- Build the constant properties for all of the species
        void buildConstProps();

        //- Record which particles are in which cell
        void buildCellOccupancy();

        //- Initialise the system
        void initialise(const IOdictionary& dsmcInitialiseDict);

        //- Calculate collisions between molecules
        void collisions();

        //- Disallow default bitwise copy construct
        DsmcCloud(const DsmcCloud&);

        //- Disallow default bitwise assignment
        void operator=(const DsmcCloud&);


public:

    // Static data members

        //- Boltzmann constant
        static scalar kb;

        //- Reference temperature for all models
        static scalar Tref;


    // Constructors

        //- Construct given name and mesh, will read Parcels from file
        DsmcCloud
        (
            const word& cloudType,
            const fvMesh& mesh
        );

        //- Construct given name, mesh, and typeIdListDict.  Used by
        // dsmcInitialse.
        DsmcCloud
        (
            const word& cloudType,
            const fvMesh& mesh,
            const IOdictionary& dsmcInitialiseDict
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

                //- Return the idList
                inline const List<word>& typeIdList() const;

                //- Return the number of real particles represented by one
                //  parcel
                inline label nParticle() const;

                //- Return the cell occupancy addressing
                inline const List<DynamicList<ParcelType*> >&
                    cellOccupancy() const;

                //- Return the sigmaTcRMax field.  non-const access to allow
                // updating.
                inline volScalarField& sigmaTcRMax();

                //- Return the collision selection remainder field.  non-const
                // access to allow updating.
                inline scalarField& collisionSelectionRemainder();

                //- Return all of the constant properties
                inline const List<typename ParcelType::constantProperties>&
                    constProps() const;

                //- Return the constant properties of the given typeId
                inline const typename ParcelType::constantProperties&
                    constProps(label typeId) const;

                //- Return refernce to the random object
                inline Random& rndGen();


            // Kinetic theory helper functions

                //- Generate a random velocity sampled from the Maxwellian speed
                // distribution
                vector equipartitionLinearVelocity
                (
                    scalar temperature,
                    scalar mass
                );

                //- From the Maxwellian distribution:

                //- Average particle speed
                inline scalar maxwellianAverageSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                //- RMS particle speed
                inline scalar maxwellianRMSSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                //- Most probable speed
                inline scalar maxwellianMostProbableSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

            //- Clear the Cloud
            inline void clear();


            // Sub-models

                //- Return reference to binary elastic collision model
                inline const BinaryCollisionModel<DsmcCloud<ParcelType> >&
                    binaryCollision() const;

                //- Return non-const reference to binary elastic collision model
                inline BinaryCollisionModel<DsmcCloud<ParcelType> >&
                    binaryCollision();

                //- Return reference to wall interaction model
                inline const WallInteractionModel<DsmcCloud<ParcelType> >&
                    wallInteraction() const;

                //- Return non-const reference to wall interaction model
                inline WallInteractionModel<DsmcCloud<ParcelType> >&
                    wallInteraction();


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
                inline const tmp<volScalarField> totalKE() const;


        // Cloud evolution functions

            //- Add new parcel
            void addNewParcel
            (
                const vector& position,
                const vector& U,
                const label cellId,
                const label typeId
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
