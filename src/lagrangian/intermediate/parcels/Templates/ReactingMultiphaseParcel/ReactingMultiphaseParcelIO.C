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

\*---------------------------------------------------------------------------*/

#include "ReactingMultiphaseParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const Cloud<ParcelType>& cloud,
    Istream& is,
    bool readFields
)
:
    ReactingParcel<ParcelType>(cloud, is, readFields),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0)
{
    if (readFields)
    {
        const ReactingCloud<ParcelType>& cR =
            dynamic_cast<const ReactingCloud<ParcelType>& >(cloud);

        const label nGas = cR.composition().gasNames().size();
        const label nLiquid = cR.composition().liquidNames().size();
        const label nSolid = cR.composition().solidNames().size();

        YGas_.setSize(nGas);
        YLiquid_.setSize(nLiquid);
        YSolid_.setSize(nSolid);

        const scalarField& YMix = this->YMixture_;
        if (is.format() == IOstream::ASCII)
        {
            is >> YGas_ >> YLiquid_ >> YSolid_;
            YGas_ /= YMix[0] + VSMALL;
            YLiquid_ /= YMix[1] + VSMALL;
            YSolid_ /= YMix[2] + VSMALL;
        }
        else
        {
            is >> YGas_ >> YLiquid_ >> YSolid_;
            YGas_ /= YMix[0] + VSMALL;
            YLiquid_ /= YMix[1] + VSMALL;
            YSolid_ /= YMix[2] + VSMALL;
        }
    }

    // Check state of Istream
    is.check
    (
        "ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel\n"
        "(\n"
        "    const Cloud<ParcelType>&,\n"
        "    Istream&,\n"
        "    bool\n"
        ")"
    );
}


template<class ParcelType>
void Foam::ReactingMultiphaseParcel<ParcelType>::readFields
(
    ReactingCloud<ParcelType>& c
)
{
    if (!c.size())
    {
        return;
    }

    ReactingParcel<ParcelType>::readFields(c);

    // Get names and sizes for each Y...
    const wordList gasNames = c.composition().gasNames();
    const wordList liquidNames = c.composition().liquidNames();
    const wordList solidNames = c.composition().solidNames();
    const label nGas = gasNames.size();
    const label nLiquid = liquidNames.size();
    const label nSolid = solidNames.size();

    // Set storage for each Y... for each parcel
    forAllIter(typename Cloud<ParcelType>, c, iter)
    {
        ReactingMultiphaseParcel<ParcelType>& p = iter();
        p.YGas_.setSize(nGas, 0.0);
        p.YLiquid_.setSize(nLiquid, 0.0);
        p.YSolid_.setSize(nSolid, 0.0);
    }

    // Populate YGas for each parcel
    forAll(gasNames, j)
    {
        IOField<scalar> YGas
        (
            c.fieldIOobject("Y" + gasNames[j], IOobject::MUST_READ)
        );

        label i = 0;
        forAllIter(typename Cloud<ParcelType>, c, iter)
        {
            ReactingMultiphaseParcel<ParcelType>& p = iter();
            p.YGas_[j] = YGas[i++]/p.YMixture()[0];
        }
    }
    // Populate YLiquid for each parcel
    forAll(liquidNames, j)
    {
        IOField<scalar> YLiquid
        (
            c.fieldIOobject("Y" + liquidNames[j], IOobject::MUST_READ)
        );

        label i = 0;
        forAllIter(typename Cloud<ParcelType>, c, iter)
        {
            ReactingMultiphaseParcel<ParcelType>& p = iter();
            p.YLiquid_[j] = YLiquid[i++]/p.YMixture()[1];
        }
    }
    // Populate YSolid for each parcel
    forAll(solidNames, j)
    {
        IOField<scalar> YSolid
        (
            c.fieldIOobject("Y" + solidNames[j], IOobject::MUST_READ)
        );

        label i = 0;
        forAllIter(typename Cloud<ParcelType>, c, iter)
        {
            ReactingMultiphaseParcel<ParcelType>& p = iter();
            p.YSolid_[j] = YSolid[i++]/p.YMixture()[2];
        }
    }
}


template<class ParcelType>
void Foam::ReactingMultiphaseParcel<ParcelType>::writeFields
(
    const ReactingCloud<ParcelType>& c
)
{
    ReactingParcel<ParcelType>::writeFields(c);

    label np =  c.size();

    // Write the composition fractions
    if (np > 0)
    {
        const wordList& gasNames = c.composition().gasNames();
        forAll(gasNames, j)
        {
            IOField<scalar> YGas
            (
                c.fieldIOobject("Y" + gasNames[j], IOobject::NO_READ),
                np
            );

            label i = 0;
            forAllConstIter(typename Cloud<ParcelType>, c, iter)
            {
                const ReactingMultiphaseParcel<ParcelType>& p0 = iter();
                YGas[i++] = p0.YGas()[j]*p0.YMixture()[0];
            }

            YGas.write();
        }
        const wordList& liquidNames = c.composition().liquidNames();
        forAll(liquidNames, j)
        {
            IOField<scalar> YLiquid
            (
                c.fieldIOobject("Y" + liquidNames[j], IOobject::NO_READ),
                np
            );

            label i = 0;
            forAllConstIter(typename Cloud<ParcelType>, c, iter)
            {
                const ReactingMultiphaseParcel<ParcelType>& p0 = iter();
                YLiquid[i++] = p0.YLiquid()[j]*p0.YMixture()[1];
            }

            YLiquid.write();
        }
        const wordList& solidNames = c.composition().solidNames();
        forAll(solidNames, j)
        {
            IOField<scalar> YSolid
            (
                c.fieldIOobject("Y" + solidNames[j], IOobject::NO_READ),
                np
            );

            label i = 0;
            forAllConstIter(typename Cloud<ParcelType>, c, iter)
            {
                const ReactingMultiphaseParcel<ParcelType>& p0 = iter();
                YSolid[i++] = p0.YSolid()[j]*p0.YMixture()[2];
            }

            YSolid.write();
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactingMultiphaseParcel<ParcelType>& p
)
{
    scalarField YGasLoc = p.YGas()*p.YMixture()[0];
    scalarField YLiquidLoc = p.YLiquid()*p.YMixture()[1];
    scalarField YSolidLoc = p.YSolid()*p.YMixture()[2];
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ReactingParcel<ParcelType>& >(p)
            << token::SPACE << YGasLoc
            << token::SPACE << YLiquidLoc
            << token::SPACE << YSolidLoc;
    }
    else
    {
        os  << static_cast<const ReactingParcel<ParcelType>& >(p);
        os << YGasLoc << YLiquidLoc << YSolidLoc;
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<\n"
        "(\n"
        "    Ostream&,\n"
        "    const ReactingMultiphaseParcel<ParcelType>&\n"
        ")"
    );

    return os;
}


// ************************************************************************* //
