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

#include "evaporationProperties.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evaporationProperties::evaporationProperties(Istream& is)
:
    name_("unknown"),
    Dab_(0.0),
    TvsPSat_(NULL)
{
    is.check
    (
        "Foam::evaporationProperties::evaporationProperties(Istream& is)"
    );

    dictionaryEntry evapInfo(dictionary::null, is);

    if (!evapInfo.isDict())
    {
        FatalErrorIn
        (
            "Foam::evaporationProperties::evaporationProperties(Istream& is)"
        )   << "Evaporation properties should be given in dictionary entries"
            << nl << exit(FatalError);
    }

    name_ = evapInfo.keyword();

    evapInfo.lookup("Dab") >> Dab_;
    TvsPSat_.reset(DataEntry<scalar>::New("TvsPSat", evapInfo.dict()).ptr());
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, evaporationProperties& ep)
{
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "("
            "Foam::Istream&,"
            "Foam::evaporationProperties&"
        ")"
    );

    dictionaryEntry evapInfo(dictionary::null, is);

    if (!evapInfo.isDict())
    {
        FatalErrorIn
        (
            "Foam::Istream& Foam::operator>>"
            "("
                "Istream& is,"
                "evaporationProperties& pp"
            ")"
        )   << "Evaporation properties should be given in dictionary entries"
            << nl << exit(FatalError);
    }

    ep.name_ = evapInfo.keyword();

    evapInfo.lookup("Dab") >> ep.Dab_;
    ep.TvsPSat_.reset(DataEntry<scalar>::New("TvsPSat", evapInfo.dict()).ptr());

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const evaporationProperties& ep)
{
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "("
            "Foam::Ostream&,"
            "const Foam::evaporationProperties&"
        ")"
    );

    os  << ep.name_ << nl << token::BEGIN_BLOCK << nl
        << incrIndent;

    os.writeKeyword("Dab") << ep.Dab_ << token::END_STATEMENT << nl;
//    os  <<  ep.TvsPSat_() << nl;

    os  << decrIndent << token::END_BLOCK << nl;

    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "("
            "Foam::Ostream&,"
            "const Foam::evaporationProperties&"
        ")"
    );

    return os;
}


// ************************************************************************* //
