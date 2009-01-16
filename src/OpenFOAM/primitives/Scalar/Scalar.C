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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const char* const Foam::pTraits<Foam::Scalar>::typeName = "scalar";
const Foam::Scalar Foam::pTraits<Foam::Scalar>::zero(0.0);
const Foam::Scalar Foam::pTraits<Foam::Scalar>::one(1.0);

const char* Foam::pTraits<Foam::Scalar>::componentNames[] = { "x" };

Foam::pTraits<Foam::Scalar>::pTraits(Istream& is)
{
    is >> p_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const Scalar s)
{
    std::ostringstream buf;
    buf << s;
    return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Scalar Foam::readScalar(Istream& is)
{
    Scalar val;
    is >> val;

    return val;
}


Foam::Istream& Foam::operator>>(Istream& is, Scalar& s)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isNumber())
    {
        s = t.number();
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, Scalar&)", is)
            << "wrong token type - expected Scalar found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, Scalar&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const Scalar s)
{
    os.write(s);
    os.check("Ostream& operator<<(Ostream&, const Scalar&)");
    return os;
}


// ************************************************************************* //
