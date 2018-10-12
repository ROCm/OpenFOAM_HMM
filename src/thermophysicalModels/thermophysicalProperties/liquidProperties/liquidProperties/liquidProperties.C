/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "liquidProperties.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquidProperties, 0);
    defineRunTimeSelectionTable(liquidProperties,);
    defineRunTimeSelectionTable(liquidProperties, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidProperties::liquidProperties
(
    scalar W,
    scalar Tc,
    scalar Pc,
    scalar Vc,
    scalar Zc,
    scalar Tt,
    scalar Pt,
    scalar Tb,
    scalar dipm,
    scalar omega,
    scalar delta
)
:
    thermophysicalProperties(W),
    Tc_(Tc),
    Pc_(Pc),
    Vc_(Vc),
    Zc_(Zc),
    Tt_(Tt),
    Pt_(Pt),
    Tb_(Tb),
    dipm_(dipm),
    omega_(omega),
    delta_(delta)
{}


Foam::liquidProperties::liquidProperties(const dictionary& dict)
:
    thermophysicalProperties(dict),
    Tc_(dict.get<label>("Tc")),
    Pc_(dict.get<label>("Pc")),
    Vc_(dict.get<label>("Vc")),
    Zc_(dict.get<label>("Zc")),
    Tt_(dict.get<label>("Tt")),
    Pt_(dict.get<label>("Pt")),
    Tb_(dict.get<label>("Tb")),
    dipm_(dict.get<label>("dipm")),
    omega_(dict.get<label>("omega")),
    delta_(dict.get<label>("delta"))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::liquidProperties> Foam::liquidProperties::New
(
    const word& name
)
{
    if (debug)
    {
        InfoInFunction << "Constructing liquidProperties" << endl;
    }

    auto cstrIter = ConstructorTablePtr_->cfind(name);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown liquidProperties type "
            << name << nl << nl
            << "Valid liquidProperties types :" << nl
            << ConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<liquidProperties>(cstrIter()());
}


Foam::autoPtr<Foam::liquidProperties> Foam::liquidProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing liquidProperties" << endl;
    }

    const word liquidType(dict.dictName());

    if (dict.found("defaultCoeffs"))
    {
        // Backward-compatibility

        if (dict.get<bool>("defaultCoeffs"))
        {
            return New(liquidType);
        }

        auto cstrIter = dictionaryConstructorTablePtr_->cfind(liquidType);

        if (!cstrIter.found())
        {
            FatalErrorInFunction
                << "Unknown liquidProperties type "
                << liquidType << nl << nl
                << "Valid liquidProperties types :" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<liquidProperties>
        (
            cstrIter()
            (
                dict.optionalSubDict(liquidType + "Coeffs")
            )
        );
    }

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(liquidType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown liquidProperties type "
            << liquidType << nl << nl
            << "Valid liquidProperties types :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<liquidProperties>(cstrIter()(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::liquidProperties::S(scalar p, scalar T) const
{
    NotImplemented;
    return 0;
}


Foam::scalar Foam::liquidProperties::pvInvert(scalar p) const
{
    // Check for critical and solid phase conditions
    if (p >= Pc_)
    {
        return Tc_;
    }
    else if (p < Pt_)
    {
        if (debug)
        {
            WarningInFunction
                << "Pressure below triple point pressure: "
                << "p = " << p << " < Pt = " << Pt_ <<  nl << endl;
        }
        return -1;
    }

    // Set initial upper and lower bounds
    scalar Thi = Tc_;
    scalar Tlo = Tt_;

    // Initialise T as boiling temperature under normal conditions
    scalar T = Tb_;

    while ((Thi - Tlo) > 1.0e-4)
    {
        if ((pv(p, T) - p) <= 0)
        {
            Tlo = T;
        }
        else
        {
            Thi = T;
        }

        T = (Thi + Tlo)*0.5;
    }

    return T;
}


void Foam::liquidProperties::readIfPresent(const dictionary &dict)
{
    thermophysicalProperties::readIfPresent(dict);
    dict.readIfPresent("Tc", Tc_);
    dict.readIfPresent("Pc", Pc_);
    dict.readIfPresent("Vc", Vc_);
    dict.readIfPresent("Zc", Zc_);
    dict.readIfPresent("Tt", Tt_);
    dict.readIfPresent("Pt", Pt_);
    dict.readIfPresent("Tb", Tb_);
    dict.readIfPresent("dipm", dipm_);
    dict.readIfPresent("omega", omega_);
    dict.readIfPresent("delta", delta_);
}


void Foam::liquidProperties::writeData(Ostream& os) const
{
    thermophysicalProperties::writeData(os);
    os  << token::SPACE
        << Tc_ << token::SPACE
        << Pc_ << token::SPACE
        << Vc_ << token::SPACE
        << Zc_ << token::SPACE
        << Tt_ << token::SPACE
        << Pt_ << token::SPACE
        << Tb_ << token::SPACE
        << dipm_ << token::SPACE
        << omega_<< token::SPACE
        << delta_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const liquidProperties& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
