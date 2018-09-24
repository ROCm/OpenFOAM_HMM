/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "cylindricalCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Issue warning if 'degrees' keyword was specified and true.
// Compatibility change after 1806.

static inline void warnCompatDegrees(const Foam::dictionary& dict)
{
    if (Pstream::parRun() ? Pstream::master() : true)
    {
        std::cerr
            << "--> FOAM IOWarning :" << nl
            << "    Found [v1806] 'degrees' keyword in dictionary \""
            << dict.name().c_str() << "\"    Ignored, now radians only." << nl
            << std::endl;
    }
}


//- Convert from Cartesian (to Cylindrical)
static inline vector fromCartesian(const vector& v)
{
    return vector(hypot(v.x(), v.y()), atan2(v.y(), v.x()), v.z());
}


//- Convert to Cartesian (from Cylindrical)
static inline vector toCartesian(const vector& v)
{
    return vector(v.x()*cos(v.y()), v.x()*sin(v.y()), v.z());
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylindricalCS::cylindricalCS()
:
    coordinateSystem()
{}


Foam::cylindricalCS::cylindricalCS
(
    const coordinateSystem& cs
)
:
    coordinateSystem(cs)
{}


Foam::cylindricalCS::cylindricalCS
(
    const word& name,
    const coordinateSystem& cs
)
:
    coordinateSystem(name, cs)
{}


Foam::cylindricalCS::cylindricalCS
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr
)
:
    coordinateSystem(name, origin, cr)
{}


Foam::cylindricalCS::cylindricalCS
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    coordinateSystem(name, origin, axis, dirn)
{}


Foam::cylindricalCS::cylindricalCS
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict)
{
    if (dict.lookupOrDefault("degrees", false))
    {
        warnCompatDegrees(dict);
    }
}


Foam::cylindricalCS::cylindricalCS
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    coordinateSystem(obr, dict)
{
    if (dict.lookupOrDefault("degrees", false))
    {
        warnCompatDegrees(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::cylindricalCS::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    return coordinateSystem::localToGlobal
    (
        toCartesian(local),
        translate
    );
}


Foam::tmp<Foam::vectorField> Foam::cylindricalCS::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    const label len = local.size();

    auto tresult = tmp<vectorField>::New(len);
    auto& result = tresult.ref();

    for (label i=0; i<len; ++i)
    {
        result[i] =
            coordinateSystem::localToGlobal
            (
                toCartesian(local[i]),
                translate
            );
    }

    return tresult;
}


Foam::vector Foam::cylindricalCS::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    return fromCartesian
    (
        coordinateSystem::globalToLocal(global, translate)
    );
}


Foam::tmp<Foam::vectorField> Foam::cylindricalCS::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    const label len = global.size();

    tmp<vectorField> tresult
    (
        coordinateSystem::globalToLocal(global, translate)
    );
    auto& result = tresult.ref();

    for (label i=0; i<len; ++i)
    {
        result[i] = fromCartesian(result[i]);
    }

    return tresult;
}


// ************************************************************************* //
