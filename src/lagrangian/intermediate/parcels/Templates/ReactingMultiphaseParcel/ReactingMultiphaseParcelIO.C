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
        const ReactingMultiphaseCloud<ParcelType>& cR =
            dynamic_cast<const ReactingMultiphaseCloud<ParcelType>& >(cloud);

        const label idGas = cR.composition().idGas();
        const label nGas = cR.composition().componentNames(idGas).size();
        const label idLiquid = cR.composition().idLiquid();
        const label nLiquid = cR.composition().componentNames(idLiquid).size();
        const label idSolid = cR.composition().idGas();
        const label nSolid = cR.composition().componentNames(idSolid).size();

        YGas_.setSize(nGas);
        YLiquid_.setSize(nLiquid);
        YSolid_.setSize(nSolid);

        if (is.format() == IOstream::ASCII)
        {
            is >> YGas_ >> YLiquid_ >> YSolid_;
        }
        else
        {
            is >> YGas_ >> YLiquid_ >> YSolid_;
        }

        // scale the mass fractions
        const scalarField& YMix = this->YMixture_;
        YGas_ /= YMix[0] + ROOTVSMALL;
        YLiquid_ /= YMix[1] + ROOTVSMALL;
        YSolid_ /= YMix[2] + ROOTVSMALL;
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
    ReactingMultiphaseCloud<ParcelType>& c
)
{
    if (!c.size())
    {
        return;
    }

    ReactingParcel<ParcelType>::readFields(c);

    // Get names and sizes for each Y...
    const label idGas = c.composition().idGas();
    const wordList gasNames = c.composition().componentNames(idGas);
    const label idLiquid = c.composition().idLiquid();
    const wordList liquidNames = c.composition().componentNames(idLiquid);
    const label idSolid = c.composition().idSolid();
    const wordList solidNames = c.composition().componentNames(idSolid);

    // Set storage for each Y... for each parcel
    forAllIter(typename Cloud<ParcelType>, c, iter)
    {
        ReactingMultiphaseParcel<ParcelType>& p = iter();
        p.YGas_.setSize(gasNames.size(), 0.0);
        p.YLiquid_.setSize(liquidNames.size(), 0.0);
        p.YSolid_.setSize(solidNames.size(), 0.0);
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
            p.YGas_[j] = YGas[i++]/(p.YMixture()[0] + ROOTVSMALL);
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
            p.YLiquid_[j] = YLiquid[i++]/(p.YMixture()[1] + ROOTVSMALL);
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
            p.YSolid_[j] = YSolid[i++]/(p.YMixture()[2] + ROOTVSMALL);
        }
    }
}


template<class ParcelType>
void Foam::ReactingMultiphaseParcel<ParcelType>::writeFields
(
    const ReactingMultiphaseCloud<ParcelType>& c
)
{
    ReactingParcel<ParcelType>::writeFields(c);

    label np =  c.size();

    // Write the composition fractions
    if (np > 0)
    {
        const label idGas = c.composition().idGas();
        const wordList gasNames = c.composition().componentNames(idGas);
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
        const label idLiquid = c.composition().idLiquid();
        const wordList liquidNames = c.composition().componentNames(idLiquid);
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
        const label idSolid = c.composition().idSolid();
        const wordList solidNames = c.composition().componentNames(idSolid);
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
        os  << static_cast<const ReactingMultiphaseParcel<ParcelType>& >(p)
            << token::SPACE << YGasLoc
            << token::SPACE << YLiquidLoc
            << token::SPACE << YSolidLoc;
    }
    else
    {
        os  << static_cast<const ReactingMultiphaseParcel<ParcelType>& >(p);
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
