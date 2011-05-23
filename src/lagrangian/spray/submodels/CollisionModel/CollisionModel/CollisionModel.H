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
    Templated collision model class

SourceFiles
    CollisionModel.C
    CollisionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef CollisionModel_H
#define CollisionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "SubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class CollisionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CollisionModel
:
    public SubModelBase<CloudType>
{
public:

    //- Runtime type information
    TypeName("collisionModel");

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

        //- Construct copy
        CollisionModel(const CollisionModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CollisionModel<CloudType> > clone() const
        {
            return autoPtr<CollisionModel<CloudType> >
            (
                new CollisionModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CollisionModel();


    //- Selector
    static autoPtr<CollisionModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        virtual bool update
        (
            const scalar dt,
            cachedRandom& rndGen,
            vector& pos1,
            scalar& m1,
            scalar& d1,
            scalar& N1,
            vector& U1,
            scalar& rho1,
            scalar& T1,
            scalarField& Y1,
            const scalar sigma1,
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
            const scalar sigma2,
            const label cellj,
            const scalar volj
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCollisionModel(CloudType)                                         \
                                                                              \
    typedef CloudType::sprayCloudType sprayCloudType;                         \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        CollisionModel<sprayCloudType>,                                       \
        0                                                                     \
    );                                                                        \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        CollisionModel<sprayCloudType>,                                       \
        dictionary                                                            \
    );


#define makeCollisionModelType(SS, CloudType)                                 \
                                                                              \
    typedef CloudType::sprayCloudType sprayCloudType;                         \
    defineNamedTemplateTypeNameAndDebug(SS<sprayCloudType>, 0);               \
                                                                              \
    CollisionModel<reactingCloudType>::                                       \
        adddictionaryConstructorToTable<SS<sprayCloudType> >                  \
            add##SS##CloudType##sprayCloudType##ConstructorToTable_;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CollisionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
