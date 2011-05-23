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
    Foam::AtomizationModel

Description
    Templated phase change model class

SourceFiles
    AtomizationModel.C
    NewAtomizationModel.C

\*---------------------------------------------------------------------------*/

#ifndef AtomizationModel_H
#define AtomizationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class AtomizationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class AtomizationModel
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
    TypeName("AtomizationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        AtomizationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        AtomizationModel(CloudType& owner);

        //- Construct from dictionary
        AtomizationModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~AtomizationModel();


    //- Selector
    static autoPtr<AtomizationModel<CloudType> > New
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

        //- initial value of liquidCore
        virtual scalar initLiquidCore() const = 0;

        //- Average temperature calculation
        scalar Taverage(const scalar& Tliq, const scalar& Tc) const;

        //- flag to indicate if chi needs to be calculated
        virtual bool calcChi() const = 0;

        virtual void update
        (
            const scalar& dt,
            scalar& d,
            scalar& liquidCore,
            scalar& tc,
            const scalar& rho,
            const scalar& mu,
            const scalar& sigma,
            const scalar& massflowRate,
            const scalar& rhoAv,
            const scalar& Urel,
            const vector& pos,
            const vector& injectionPos,
            const scalar& pAmbient,
            const scalar& chi,
            Random& rndGen
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeAtomizationModel(CloudType)                                       \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(AtomizationModel<CloudType>, 0);      \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        AtomizationModel<CloudType>,                                          \
        dictionary                                                            \
    );

#define makeAtomizationModelThermoType(SS, CloudType, ParcelType, ThermoType) \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        SS<CloudType<ParcelType<ThermoType> > >,                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    AtomizationModel<CloudType<ParcelType<ThermoType> > >::                   \
        adddictionaryConstructorToTable                                       \
            <SS<CloudType<ParcelType<ThermoType> > > >                        \
            add##SS##CloudType##ParcelType##ThermoType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "AtomizationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
