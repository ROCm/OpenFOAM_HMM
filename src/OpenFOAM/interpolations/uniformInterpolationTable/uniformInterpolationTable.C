/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "uniformInterpolationTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uniformInterpolationTable, 0);
}

// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

void Foam::uniformInterpolationTable::checkTable() const
{
    if (size() < 2)
    {
        FatalErrorIn("uniformInterpolationTable::checkTable()")
            << "Table " << name() << ": must have at least 2 values." << nl
            << "Table size = " << size() << nl
            << "    min, interval width = " << x0_ << ", " << dx_ << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformInterpolationTable::uniformInterpolationTable
(
    const IOobject& io,
    bool readFields
)
:
    IOobject(io),
    List<scalar>(2, 0.0),
    x0_(0.0),
    dx_(1.0),
    log10_(false),
    bound_(false)
{
    if (readFields)
    {
        IOdictionary dict(io);

        dict.lookup("data") >> *this;
        dict.lookup("x0") >> x0_;
        dict.lookup("dx") >> dx_;
        dict.lookup("log10") >> log10_;
        dict.lookup("bound") >> bound_;
    }

    checkTable();
}


Foam::uniformInterpolationTable::uniformInterpolationTable
(
    const word& tableName,
    const objectRegistry& db,
    const dictionary& dict,
    const bool initialiseOnly
)
:
    IOobject
    (
        tableName,
        db.time().constant(),
        db,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false // if used in BCs, could be used by multiple patches
    ),
    List<scalar>(2, 0.0),
    x0_(readScalar(dict.lookup("x0"))),
    dx_(readScalar(dict.lookup("dx"))),
    log10_(dict.lookup("log10")),
    bound_(dict.lookup("bound"))
{
    if (initialiseOnly)
    {
        scalar xMax = readScalar(dict.lookup("xMax"));
        label nIntervals = static_cast<label>(xMax - x0_)/dx_ + 1;
        setSize(nIntervals);
    }
    else
    {
        dict.lookup("data") >> *this;
    }

    checkTable();
}


Foam::uniformInterpolationTable::uniformInterpolationTable(const uniformInterpolationTable& uit)
:
    IOobject(uit),
    List<scalar>(uit),
    x0_(uit.x0_),
    dx_(uit.dx_),
    log10_(uit.log10_),
    bound_(uit.bound_)
{
    checkTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uniformInterpolationTable::~uniformInterpolationTable()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::uniformInterpolationTable::interpolate(scalar x) const
{
    if (bound_)
    {
        x = max(min(xMax() - SMALL*dx_, x), x0_);
    }
    else
    {
        if (x < x0_)
        {
            FatalErrorIn("uniformInterpolationTable::interpolate(scalar x)")
                << "Supplied value is less than minimum table value:" << nl
                << "xMin=" << x0_ << ", xMax=" << xMax() << ", x=" << x << nl
                << exit(FatalError);
        }

        if (x > xMax())
        {
            FatalErrorIn("uniformInterpolationTable::interpolate(scalar x)")
                << "Supplied value is greater than maximum table value:" << nl
                << "xMin=" << x0_ << ", xMax=" << xMax() << ", x=" << x << nl
                << exit(FatalError);
        }
    }

    label i = static_cast<label>((x - x0_)/dx_);

    scalar xLo = x0_ + i*dx_;

    scalar fx = (x - xLo)/dx_*(operator[](i+1) - operator[](i)) + operator[](i);

    if (debug)
    {
        Info<< "Table: " << name() << ", x=" << x
            << ", x_lo=" << xLo << ", x_hi=" << xLo + dx_
            << ", f(x_lo)=" << operator[](i) << ", f(x_hi)=" << operator[](i+1)
            << ", f(x)=" << fx << endl;
    }

    return fx;
}


Foam::scalar Foam::uniformInterpolationTable::interpolateLog10(scalar x) const
{
    if (log10_)
    {
        if (x > 0)
        {
            x = ::log10(x);
        }
        else if (bound_ && (x <= 0))
        {
            x = x0_;
        }
        else
        {
            FatalErrorIn
            (
                "uniformInterpolationTable::interpolateLog10(scalar x)"
            )   << "Table " << name() << nl
                << "Supplied value must be greater than 0 when in log10 mode"
                << nl << "x=" << x << nl << exit(FatalError);
        }
    }

    return interpolate(x);
}


void Foam::uniformInterpolationTable::write() const
{
    IOdictionary dict(*this);

    dict.add("data", static_cast<const List<scalar>&>(*this));
    dict.add("x0", x0_);
    dict.add("dx", dx_);
    dict.add("log10", log10_);
    dict.add("bound", bound_);

    dict.regIOobject::write();
}


// ************************************************************************* //
