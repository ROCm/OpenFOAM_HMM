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

#include "fileName.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fileName::fileName(Istream& is)
:
    string()
{
    is >> *this;
}


Foam::Istream& Foam::operator>>(Istream& is, fileName& fn)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isString())
    {
        fn = t.stringToken();
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, fileName&)", is)
            << "wrong token type - expected string, found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    fn.stripInvalid();

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, fileName&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const fileName& fn)
{
    os.write(fn);
    os.check("Ostream& operator<<(Ostream&, const fileName&)");
    return os;
}


// ************************************************************************* //


