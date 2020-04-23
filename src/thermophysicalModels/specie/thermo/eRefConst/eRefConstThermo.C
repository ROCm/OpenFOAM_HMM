/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "eRefConstThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::eRefConstThermo<EquationOfState>::eRefConstThermo(const dictionary& dict)
:
    EquationOfState(dict),
    Cv_(readCoeff("Cv", dict)),
    Hf_(readCoeff("Hf", dict)),
    Tref_(readCoeff("Tref", dict)),
    Eref_(readCoeff("Eref", dict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::eRefConstThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    {
        os.beginBlock("thermodynamics");
        os.writeEntry("Cv", Cv_);
        os.writeEntry("Hf", Hf_);
        os.writeEntry("Tref", Tref_);
        os.writeEntry("Eref", Eref_);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const eRefConstThermo<EquationOfState>& ct
)
{
    ct.write(os);
    return os;
}


// ************************************************************************* //
