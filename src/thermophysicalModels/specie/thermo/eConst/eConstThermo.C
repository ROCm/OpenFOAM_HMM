/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

#include "eConstThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::eConstThermo<EquationOfState>::eConstThermo(Istream& is)
:
    EquationOfState(is),
    Cv_(readScalar(is)),
    Hf_(readScalar(is))
{
    is.check("eConstThermo<EquationOfState>::eConstThermo(Istream&)");
}


template<class EquationOfState>
Foam::eConstThermo<EquationOfState>::eConstThermo(const dictionary& dict)
:
    EquationOfState(dict),
    Cv_(readScalar(dict.lookup("Cv"))),
    Hf_(readScalar(dict.lookup("Hf")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::eConstThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);
    os.writeKeyword("Cv") << Cv_ << token::END_STATEMENT << nl;
    os.writeKeyword("Hf") << Hf_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const eConstThermo<EquationOfState>& ct
)
{
    os  << static_cast<const EquationOfState&>(ct) << tab
        << ct.Cv_ << tab << ct.Hf_;

    os.check("Ostream& operator<<(Ostream&, const eConstThermo&)");
    return os;
}


// ************************************************************************* //
