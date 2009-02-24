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
    Foam::BinaryElasticCollisionModel

Description
    Templated DSMC particle collision class

SourceFiles
    BinaryElasticCollisionModel.C
    NewBinaryElasticCollisionModel.C

\*---------------------------------------------------------------------------*/

#ifndef BinaryElasticCollisionModel_H
#define BinaryElasticCollisionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BinaryElasticCollisionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BinaryElasticCollisionModel
{

    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        // reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;


public:

    //- Runtime type information
    TypeName("BinaryElasticCollisionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        BinaryElasticCollisionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct from components
        BinaryElasticCollisionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    // Destructor

        virtual ~BinaryElasticCollisionModel();


    // Selector

        static autoPtr<BinaryElasticCollisionModel<CloudType> > New
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Access

        //- Return the owner cloud object
        const CloudType& owner() const;

        //- Return the dictionary
        const dictionary& dict() const;

        //- Return the coefficients dictionary
        const dictionary& coeffDict() const;


    // Member Functions

        //- Apply collision
        virtual void collide
        (
            label idA,
            label idB,
            vector& UA,
            vector& UB
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeBinaryElasticCollisionModel(CloudType)                            \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        BinaryElasticCollisionModel<CloudType>,                               \
        0                                                                     \
    );                                                                        \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        BinaryElasticCollisionModel<CloudType>,                               \
        dictionary                                                            \
    );


#define makeBinaryElasticCollisionModelType(SS, CloudType, ParcelType)        \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType<ParcelType> >, 0);       \
                                                                              \
    BinaryElasticCollisionModel<CloudType<ParcelType> >::                     \
        adddictionaryConstructorToTable<SS<CloudType<ParcelType> > >          \
            add##SS##CloudType##ParcelType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BinaryElasticCollisionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
