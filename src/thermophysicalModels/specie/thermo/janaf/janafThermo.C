/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "janafThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class equationOfState>
void Foam::janafThermo<equationOfState>::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorIn("janafThermo<equationOfState>::check()")
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalIOError);
    }

    if (Tcommon_ <= Tlow_)
    {
        FatalErrorIn("janafThermo<equationOfState>::check()")
            << "Tcommon(" << Tcommon_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalIOError);
    }

    if (Tcommon_ > Thigh_)
    {
        FatalErrorIn("janafThermo<equationOfState>::check()")
            << "Tcommon(" << Tcommon_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
Foam::janafThermo<equationOfState>::janafThermo(Istream& is)
:
    equationOfState(is),
    Tlow_(readScalar(is)),
    Thigh_(readScalar(is)),
    Tcommon_(readScalar(is))
{
    checkInputData();

    forAll(highCpCoeffs_, i)
    {
        is >> highCpCoeffs_[i];
    }

    forAll(lowCpCoeffs_, i)
    {
        is >> lowCpCoeffs_[i];
    }

    // Check state of Istream
    is.check("janafThermo::janafThermo(Istream& is)");
}


template<class equationOfState>
Foam::janafThermo<equationOfState>::janafThermo(const dictionary& dict)
:
    equationOfState(dict),
    Tlow_(readScalar(dict.lookup("Tlow"))),
    Thigh_(readScalar(dict.lookup("Thigh"))),
    Tcommon_(readScalar(dict.lookup("Tcommon"))),
    highCpCoeffs_(dict.lookup("highCpCoeffs")),
    lowCpCoeffs_(dict.lookup("lowCpCoeffs"))
{
    checkInputData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState>
void Foam::janafThermo<equationOfState>::write(Ostream& os) const
{
    equationOfState::write(os);
    os.writeKeyword("Tlow") << Tlow_ << token::END_STATEMENT << endl;
    os.writeKeyword("Thigh") << Thigh_ << token::END_STATEMENT << endl;
    os.writeKeyword("Tcommon") << Tcommon_ << token::END_STATEMENT << endl;
    os.writeKeyword("highCpCoeffs") << highCpCoeffs_ << token::END_STATEMENT
        << endl;
    os.writeKeyword("lowCpCoeffs") << lowCpCoeffs_ << token::END_STATEMENT
        << endl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class equationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const janafThermo<equationOfState>& jt
)
{
    os  << static_cast<const equationOfState&>(jt) << nl
        << "    " << jt.Tlow_
        << tab << jt.Thigh_
        << tab << jt.Tcommon_;

    os << nl << "    ";

    forAll(jt.highCpCoeffs_, i)
    {
        os << jt.highCpCoeffs_[i] << ' ';
    }

    os << nl << "    ";

    forAll(jt.lowCpCoeffs_, i)
    {
        os << jt.lowCpCoeffs_[i] << ' ';
    }

    os << endl;

    os.check
    (
        "operator<<(Ostream& os, const janafThermo<equationOfState>& jt)"
    );

    return os;
}


// ************************************************************************* //
