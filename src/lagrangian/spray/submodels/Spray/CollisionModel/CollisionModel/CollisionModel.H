/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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
    Foam::CollisionModel

Description
    Templated phase change model class

SourceFiles
    CollisionModel.C
    NewCollisionModel.C

\*---------------------------------------------------------------------------*/

#ifndef CollisionModel_H
#define CollisionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class CollisionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CollisionModel
{

protected:

    // Protected data

        //- The cloud dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficient dictionary
        const dictionary coeffDict_;

public:

    //- Runtime type information
    TypeName("CollisionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CollisionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        CollisionModel(CloudType& owner);

        //- Construct from dictionary
        CollisionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~CollisionModel();


    //- Selector
    static autoPtr<CollisionModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Access

        //- Return the owner cloud object
        const CloudType& owner() const;

        //- Return the cloud dictionary
        const dictionary& dict() const;

        //- Return the coefficient dictionary
        const dictionary& coeffDict() const;


    // Member Functions

        virtual bool active() const = 0;

        virtual bool update
        (
            const scalar& dt,
            Random& rndGen,
            vector& pos1,
            scalar& m1,
            scalar& d1,
            scalar& N1,
            vector& U1,
            scalar& rho1,
            scalar& T1,
            scalarField& Y1,
            const scalar& sigma1,
            const label celli,
            const scalar voli,
            vector& pos2,
            scalar& m2,
            scalar& d2,
            scalar& N2,
            vector& U2,
            scalar& rho2,
            scalar& T2,
            scalarField& Y2,
            const scalar& sigma2,
            const label cellj,
            const scalar volj
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCollisionModel(CloudType)                                         \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(CollisionModel<CloudType>, 0);        \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        CollisionModel<CloudType>,                                            \
        dictionary                                                            \
    );

#define makeCollisionModelThermoType(SS, CloudType, ParcelType, ThermoType)   \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        SS<CloudType<ParcelType<ThermoType> > >,                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    CollisionModel<CloudType<ParcelType<ThermoType> > >::                     \
        adddictionaryConstructorToTable                                       \
            <SS<CloudType<ParcelType<ThermoType> > > >                        \
            add##SS##CloudType##ParcelType##ThermoType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CollisionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
