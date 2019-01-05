/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "ReactingHeterogeneousParcel.H"
#include "IOstreams.H"
#include "HeterogeneousReactingModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ReactingHeterogeneousParcel<ParcelType>::propertyList_ =
    Foam::ReactingHeterogeneousParcel<ParcelType>::propertyList();

template<class ParcelType>
Foam::string Foam::ReactingHeterogeneousParcel<ParcelType>::propertyTypes_ =
    Foam::ReactingHeterogeneousParcel<ParcelType>::propertyTypes();

template<class ParcelType>
const std::size_t Foam::ReactingHeterogeneousParcel<ParcelType>::sizeofFields
(
    sizeof(scalar)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingHeterogeneousParcel<ParcelType>::ReactingHeterogeneousParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    F_(0),
    canCombust_(1)
{
    if (readFields)
    {
        DynamicList<scalar> F;

        is >> F;

        F_.transfer(F);
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingHeterogeneousParcel<ParcelType>::readFields(CloudType& c)
{
    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingHeterogeneousParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    bool valid = c.size();

    ParcelType::readFields(c);

    IOField<scalar> mass0
    (
        c.fieldIOobject("mass0", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, mass0);

    label i = 0;
    for (ReactingHeterogeneousParcel<ParcelType>& p : c)
    {
        p.mass0() = mass0[i];
        ++i;
    }

    const label idSolid = compModel.idSolid();
    const wordList& solidNames = compModel.componentNames(idSolid);

    // WIP until find a way to get info from Reacting Heterogeneous model
    label nF(1);

    // Set storage for each Y... for each parcel
    for (ReactingHeterogeneousParcel<ParcelType>& p : c)
    {
        p.Y().setSize(solidNames.size(), Zero);
        p.F_.setSize(nF, Zero);
    }

    for (label i = 0; i < nF; i++)
    {
        // Read F
        IOField<scalar> F
        (
            c.fieldIOobject
            (
                "F" + name(i),
                IOobject::MUST_READ
            ),
            valid
        );

        label j = 0;
        for (ReactingHeterogeneousParcel<ParcelType>& p : c)
        {
            p.F_[i] = F[j];
            ++j;
        }
    }


    forAll(solidNames, j)
    {
        IOField<scalar> Y
        (
            c.fieldIOobject
            (
                "Y" + solidNames[j],
                IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        for (ReactingHeterogeneousParcel<ParcelType>& p : c)
        {
            p.Y()[j] = Y[i];
            ++i;
        }
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingHeterogeneousParcel<ParcelType>::writeFields
(
    const CloudType& c
)
{
    ParcelType::writeFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingHeterogeneousParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    // Skip Reacting layer. This class takes charge of
    // writing Ysolids and F
    ThermoParcel<KinematicParcel<particle>>::writeFields(c);

    label np = c.size();

    IOField<scalar> mass0(c.fieldIOobject("mass0", IOobject::NO_READ), np);

    label nF = 0;
    label i = 0;
    for (const ReactingHeterogeneousParcel<ParcelType>& p : c)
    {
        mass0[i] = p.mass0_;
        if (!i)
        {
            // Assume same size throughout
            nF = p.F().size();
        }
        ++i;
    }
    mass0.write(np > 0);

    for (label i = 0; i < nF; i++)
    {
        IOField<scalar> F
        (
            c.fieldIOobject
            (
                "F" + name(i),
                IOobject::NO_READ
            ),
            np
        );

        for (const ReactingHeterogeneousParcel<ParcelType>& p0 : c)
        {
            F = p0.F()[i];
        }

        F.write(np > 0);
    }

    const label idSolid = compModel.idSolid();
    const wordList& solidNames = compModel.componentNames(idSolid);

    forAll(solidNames, j)
    {
        IOField<scalar> Y
        (
            c.fieldIOobject
            (
                "Y" + solidNames[j],
                IOobject::NO_READ
            ),
            np
        );

        label i = 0;
        for (const ReactingHeterogeneousParcel<ParcelType>& p0 : c)
        {
            Y[i] = p0.Y()[j];
            ++i;
        }

        Y.write(np > 0);
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingHeterogeneousParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingHeterogeneousParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    const CompositionType& compModel,
    objectRegistry& obr
)
{
    //ParcelType::writeObjects(c, obr);
    // Skip Reacting layer
    ThermoParcel<KinematicParcel<particle>>::writeObjects(c, obr);

    label np = c.size();

    // WIP
    label nF = 0;
    for (const ReactingHeterogeneousParcel<ParcelType>& p0 : c)
    {
        nF = p0.F().size();
        break;
    }

    if (np > 0)
    {
        for (label i = 0; i < nF; i++)
        {
            const word fieldName = "F" + name(i);
            IOField<scalar>& F
            (
                cloud::createIOField<scalar>(fieldName, np, obr)
            );

            label j = 0;
            for (const ReactingHeterogeneousParcel<ParcelType>& p0 : c)
            {
                F[j] = p0.F()[i];
                ++j;
            }
        }

        const label idSolid = compModel.idSolid();
        const wordList& solidNames = compModel.componentNames(idSolid);
        forAll(solidNames, j)
        {
            const word fieldName = "Y" + solidNames[j];
            IOField<scalar>& Y
            (
                cloud::createIOField<scalar>(fieldName, np, obr)
            );

            label i = 0;
            for (const ReactingHeterogeneousParcel<ParcelType>& p0 : c)
            {
                Y[i] = p0.Y()[j];
                ++i;
            }
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactingHeterogeneousParcel<ParcelType>& p
)
{
    scalarField F(p.F());
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << F;
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os  << F ;
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
