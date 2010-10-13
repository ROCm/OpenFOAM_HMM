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

#include "constTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
constTransport<thermo>::constTransport(Istream& is)
:
    thermo(is),
    Mu_(readScalar(is)),
    rPr_(1.0/readScalar(is))
{
    is.check("constTransport::constTransport(Istream& is)");
}


template<class thermo>
constTransport<thermo>::constTransport(const dictionary& dict)
:
    thermo(dict),
    Mu_(readScalar(dict.lookup("Mu"))),
    rPr_(1.0/readScalar(dict.lookup("Pr")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
void constTransport<thermo>::constTransport::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;
    thermo::write(os);
    os.writeKeyword("Mu") << Mu_ << token::END_STATEMENT << nl;
    os.writeKeyword("Pr") << Mu_ << token::END_STATEMENT << nl;
    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class thermo>
Ostream& operator<<(Ostream& os, const constTransport<thermo>& ct)
{
    operator<<(os, static_cast<const thermo&>(ct));
    os << tab << ct.Mu_ << tab << 1.0/ct.rPr_;

    os.check("Ostream& operator<<(Ostream& os, const constTransport& ct)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
