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
    Foam::InjectionModel


Description
    Templated injection model class

SourceFiles
    InjectionModel.C
    NewInjectionModel.C

\*---------------------------------------------------------------------------*/

#ifndef InjectionModel_H
#define InjectionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class InjectionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InjectionModel
{

    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        // Reference to the owner cloud class
        CloudType& owner_;


public:

    //- Runtime type information
    TypeName("InjectionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        InjectionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct from dictionary
        InjectionModel
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Destructor

        virtual ~InjectionModel();


    // Selector

        static autoPtr<InjectionModel<CloudType> > New
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Access

        //- Return the owner cloud object
        const CloudType& owner() const;

        //- Return the dictionary
        const dictionary& dict() const;


    // Member Functions

        //- Flag to indicate whether model activates injection model
        virtual bool active() const = 0;

        //- Return the start-of-injection time
        virtual scalar timeStart() const = 0;

        //- Return the end-of-injection time
        virtual scalar timeEnd() const = 0;

        //- Return the injection position
        virtual vector position
        (
            const label iParcel,
            const scalar time,
            const polyMeshInfo& meshInfo,
            Random& rndGen
        ) const = 0;

        //- Return the number of parcels to intruduce between two times
        virtual label nParcelsToInject
        (
            const label nInjections_,
            const scalar time0,
            const scalar time1
        ) const = 0;

        //- Return the volume of parcels to introduce between two times
        virtual scalar volume
        (
            const scalar time0,
            const scalar time1,
            const polyMeshInfo& meshInfo
        ) const = 0;

        //- Return the volume fraction to introduce between two times
        virtual scalar volumeFraction
        (
            const scalar time0,
            const scalar time1
        ) const = 0;

        //- Return the diameter of the parcel to introduce at a time
        virtual scalar d0
        (
            const label iParcel,
            const scalar time
        ) const = 0;

        //- Return the velocity of the parcel to introduce at a time
        virtual vector velocity
        (
            const label iParcel,
            const scalar time,
            const polyMeshInfo& meshInfo,
            Random& rndGen
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeInjectionModel(CloudType)                                         \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(InjectionModel<CloudType>, 0);        \
                                                                              \
    defineTemplateRunTimeSelectionTable(InjectionModel<CloudType>, dictionary);


#define makeInjectionModelType(SS, CloudType, ParcelType)                     \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType<ParcelType> >, 0);       \
                                                                              \
    InjectionModel<CloudType<ParcelType> >::                                  \
        adddictionaryConstructorToTable<SS<CloudType<ParcelType> > >          \
            add##SS##CloudType##ParcelType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "InjectionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
