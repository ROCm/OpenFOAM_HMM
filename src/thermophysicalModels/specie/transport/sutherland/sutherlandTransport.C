/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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

#include "sutherlandTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::sutherlandTransport<Thermo>::sutherlandTransport(Istream& is)
:
    Thermo(is),
    As_(readScalar(is)),
    Ts_(readScalar(is))
{
    is.check("sutherlandTransport<Thermo>::sutherlandTransport(Istream&)");
}


template<class Thermo>
Foam::sutherlandTransport<Thermo>::sutherlandTransport(const dictionary& dict)
:
    Thermo(dict),
    As_(readScalar(dict.lookup("As"))),
    Ts_(readScalar(dict.lookup("Ts")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::sutherlandTransport<Thermo>::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;
    Thermo::write(os);
    os.writeKeyword("As") << As_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ts") << Ts_ << token::END_STATEMENT << nl;
    os  << decrIndent << token::END_BLOCK << nl;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os, const sutherlandTransport<Thermo>& st
)
{
    os << static_cast<const Thermo&>(st) << tab << st.As_ << tab << st.Ts_;

    os.check
    (
        "Ostream& operator<<(Ostream&, const sutherlandTransport<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
