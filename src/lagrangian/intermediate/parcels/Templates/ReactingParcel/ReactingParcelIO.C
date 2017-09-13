/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "ReactingParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ReactingParcel<ParcelType>::propertyList_ =
    Foam::ReactingParcel<ParcelType>::propertyList();

template<class ParcelType>
Foam::string Foam::ReactingParcel<ParcelType>::propertyTypes_ =
    Foam::ReactingParcel<ParcelType>::propertyTypes();

template<class ParcelType>
const std::size_t Foam::ReactingParcel<ParcelType>::sizeofFields
(
    sizeof(scalar)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    mass0_(0.0),
    Y_(0)
{
    if (readFields)
    {
        DynamicList<scalar> Ymix;

        if (is.format() == IOstream::ASCII)
        {
            is >> mass0_ >> Ymix;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&mass0_), sizeofFields);
            is >> Ymix;
        }

        Y_.transfer(Ymix);
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingParcel<ParcelType>::readFields(CloudType& c)
{
    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingParcel<ParcelType>::readFields
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
    forAllIter(typename Cloud<ReactingParcel<ParcelType>>, c, iter)
    {
        ReactingParcel<ParcelType>& p = iter();
        p.mass0_ = mass0[i++];
    }

    // Get names and sizes for each Y...
    const wordList& phaseTypes = compModel.phaseTypes();
    const label nPhases = phaseTypes.size();
    wordList stateLabels(nPhases, "");
    if (compModel.nPhase() == 1)
    {
        stateLabels = compModel.stateLabels()[0];
    }


    // Set storage for each Y... for each parcel
    forAllIter(typename Cloud<ReactingParcel<ParcelType>>, c, iter)
    {
        ReactingParcel<ParcelType>& p = iter();
        p.Y_.setSize(nPhases, 0.0);
    }

    // Populate Y for each parcel
    forAll(phaseTypes, j)
    {
        IOField<scalar> Y
        (
            c.fieldIOobject
            (
                "Y" + phaseTypes[j] + stateLabels[j],
                 IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        forAllIter(typename Cloud<ReactingParcel<ParcelType>>, c, iter)
        {
            ReactingParcel<ParcelType>& p = iter();
            p.Y_[j] = Y[i++];
        }
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c);

    const label np = c.size();

    {
        IOField<scalar> mass0(c.fieldIOobject("mass0", IOobject::NO_READ), np);

        label i = 0;
        forAllConstIter(typename Cloud<ReactingParcel<ParcelType>>, c, iter)
        {
            const ReactingParcel<ParcelType>& p = iter();
            mass0[i++] = p.mass0_;
        }
        mass0.write(np > 0);

        // Write the composition fractions
        const wordList& phaseTypes = compModel.phaseTypes();
        wordList stateLabels(phaseTypes.size(), "");
        if (compModel.nPhase() == 1)
        {
            stateLabels = compModel.stateLabels()[0];
        }

        forAll(phaseTypes, j)
        {
            IOField<scalar> Y
            (
                c.fieldIOobject
                (
                    "Y" + phaseTypes[j] + stateLabels[j],
                    IOobject::NO_READ
                ),
                np
            );
            label i = 0;
            forAllConstIter
            (
                typename Cloud<ReactingParcel<ParcelType>>,
                c,
                iter
            )
            {
                const ReactingParcel<ParcelType>& p = iter();
                Y[i++] = p.Y()[j];
            }

            Y.write(np > 0);
        }
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    const CompositionType& compModel,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    label np = c.size();

    if (np > 0)
    {
        IOField<scalar>& mass0(cloud::createIOField<scalar>("mass0", np, obr));

        label i = 0;
        forAllConstIter(typename Cloud<ReactingParcel<ParcelType>>, c, iter)
        {
            const ReactingParcel<ParcelType>& p = iter();
            mass0[i++] = p.mass0_;
        }

        // Write the composition fractions
        const wordList& phaseTypes = compModel.phaseTypes();
        wordList stateLabels(phaseTypes.size(), "");
        if (compModel.nPhase() == 1)
        {
            stateLabels = compModel.stateLabels()[0];
        }

        forAll(phaseTypes, j)
        {
            const word fieldName = "Y" + phaseTypes[j] + stateLabels[j];
            IOField<scalar>& Y
            (
                cloud::createIOField<scalar>(fieldName, np, obr)
            );

            label i = 0;
            forAllConstIter
            (
                typename Cloud<ReactingParcel<ParcelType>>,
                c,
                iter
            )
            {
                const ReactingParcel<ParcelType>& p = iter();
                Y[i++] = p.Y()[j];
            }
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactingParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.mass0()
            << token::SPACE << p.Y();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.mass0_),
            ReactingParcel<ParcelType>::sizeofFields
        );
        os  << p.Y();
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
