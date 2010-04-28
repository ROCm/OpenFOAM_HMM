/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "vectorTensorTransform.H"
#include "IOstreams.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::vectorTensorTransform::typeName =
"vectorTensorTransform";

const Foam::vectorTensorTransform Foam::vectorTensorTransform::zero
(
    vector::zero,
    tensor::zero
);


const Foam::vectorTensorTransform Foam::vectorTensorTransform::I
(
    vector::zero,
    sphericalTensor::I
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vectorTensorTransform::vectorTensorTransform(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const vectorTensorTransform& s)
{
    OStringStream buf;

    buf << '(' << s.t() << ',' << s.R() << ')';

    return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, vectorTensorTransform& tr)
{
    // Read beginning of vectorTensorTransform
    is.readBegin("vectorTensorTransform");

    is  >> tr.t() >> tr.R();

    // Read end of vectorTensorTransform
    is.readEnd("vectorTensorTransform");

    // Check state of Istream
    is.check("operator>>(Istream&, vectorTensorTransform&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const vectorTensorTransform& tr)
{
    os  << token::BEGIN_LIST
        << tr.t() << token::SPACE << tr.R()
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
