/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::SurfaceFilmModel

Description
    Templated wall surface film model class.

SourceFiles
    SurfaceFilmModel.C
    SurfaceFilmModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef SurfaceFilmModel_H
#define SurfaceFilmModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
namespace surfaceFilmModels
{
    class surfaceFilmModel;
}

class mapDistribute;

/*---------------------------------------------------------------------------*\
                      Class SurfaceFilmModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SurfaceFilmModel
{
protected:

    // Protected data

        //- Convenience typedef to the cloud's parcel type
        typedef typename CloudType::parcelType parcelType;

        //- The cloud dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& owner_;

        //- Gravitational acceleration constant
        const dimensionedVector& g_;

        //- The coefficients dictionary
        const dictionary coeffDict_;

        //- Active flag
        bool active_;


        // Cached injector fields per film patch

            //- Injector cell / patch face
            labelList injectorCellsPatch_;

            //- Parcel mass / patch face
            scalarList massParcelPatch_;

            //- Parcel diameter / patch face
            scalarList diameterParcelPatch_;

            //- Film velocity / patch face
            List<vector> UFilmPatch_;

            //- Film density / patch face
            scalarList rhoFilmPatch_;


        // Counters

            //- Number of parcels transferred to the film model
            label nParcelsTransferred_;

            //- Number of parcels injected from the film model
            label nParcelsInjected_;


    // Protected functions

        //- Cache the film fields in preparation for injection
        virtual void cacheFilmFields
        (
            const label filmPatchI,
            const mapDistribute& distMap,
            const surfaceFilmModels::surfaceFilmModel& filmModel
        );

        //- Set the individual parcel properties
        virtual void setParcelProperties
        (
            parcelType& p,
            const label filmFaceI
        ) const;


public:

    //- Runtime type information
    TypeName("SurfaceFilmModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        SurfaceFilmModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner,
            const dimensionedVector& g
        ),
        (dict, owner, g)
    );


    // Constructors

        //- Construct null from owner
        SurfaceFilmModel(CloudType& owner);

        //- Construct from dictionary
        SurfaceFilmModel
        (
            const dictionary& dict,
            CloudType& owner,
            const dimensionedVector& g,
            const word& type
        );


    //- Destructor
    virtual ~SurfaceFilmModel();


    //- Selector
    static autoPtr<SurfaceFilmModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner,
        const dimensionedVector& g
    );


    // Member Functions

        // Access

            //- Return the owner cloud dictionary
            inline const dictionary& dict() const;

            //- Return const access the owner cloud object
            inline const CloudType& owner() const;

            //- Return non-const access the owner cloud object for manipulation
            inline CloudType& owner();

            //- Return gravitational acceleration constant
            inline const dimensionedVector& g() const;

            //- Return the coefficients dictionary
            inline const dictionary& coeffDict() const;

            //- Return const access to the number of parcels transferred to the
            //  film model
            inline label nParcelsTransferred() const;

            //- Return non-const access to the number of parcels transferred to
            //  the film model
            inline label& nParcelsTransferred();

            //- Return const access to the number of parcels injected from the
            //  film model
            inline label nParcelsInjected() const;

            //- Return non-const access to the number of parcels injected from
            //  the film model
            inline label& nParcelsInjected();


        // Member Functions

            //- Flag to indicate whether model activates the surface film model
            virtual bool active() const = 0;

            //- Transfer parcel from cloud to surface film
            //  Returns true if parcel is to be transferred
            virtual bool transferParcel
            (
                const parcelType& p,
                const label patchI
            ) = 0;

            //- Inject parcels into the cloud
            template<class TrackData>
            void inject(TrackData& td);


        // I-O

            //- Write surface film info to stream
            virtual void info(Ostream& os) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSurfaceFilmModel(CloudType)                                       \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SurfaceFilmModel<CloudType>, 0);      \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        SurfaceFilmModel<CloudType>,                                          \
        dictionary                                                            \
    );


#define makeSurfaceFilmModelType(SS, CloudType, ParcelType)                   \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType<ParcelType> >, 0);       \
                                                                              \
    SurfaceFilmModel<CloudType<ParcelType> >::                                \
        adddictionaryConstructorToTable<SS<CloudType<ParcelType> > >          \
            add##SS##CloudType##ParcelType##ConstructorToTable_;


#define makeSurfaceFilmModelThermoType(SS, CloudType, ParcelType, ThermoType) \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        SS<CloudType<ParcelType<ThermoType> > >,                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    SurfaceFilmModel<CloudType<ParcelType<ThermoType> > >::                   \
        adddictionaryConstructorToTable                                       \
        <SS<CloudType<ParcelType<ThermoType> > > >                            \
            add##SS##CloudType##ParcelType##ThermoType##ConstructorToTable_;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SurfaceFilmModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SurfaceFilmModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
