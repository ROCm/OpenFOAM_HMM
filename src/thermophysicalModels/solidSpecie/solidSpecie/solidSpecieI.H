/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "solidSpecie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline solidSpecie::solidSpecie
(
    const word& name,
    const scalar nMoles,
    const scalar molWeight
)
:
    specie(name, nMoles, scalar(1))
{}


inline solidSpecie::solidSpecie
(
    const scalar nMoles, const scalar molWeight
)
:
    specie(nMoles, molWeight)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline solidSpecie::solidSpecie(const solidSpecie& st)
:
    specie(st.name(), st.nMoles(), scalar(1))
{}


inline solidSpecie::solidSpecie(const word& name, const solidSpecie& st)
:
    specie(name, st.nMoles(), scalar(1))
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void solidSpecie::operator=(const solidSpecie& st)
{
    specie::operator=(st);
}


inline void solidSpecie::operator+=(const solidSpecie& st)
{
    specie::operator+=(st);
}


inline void solidSpecie::operator-=(const solidSpecie& st)
{
    specie::operator-=(st);
}


inline void solidSpecie::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline solidSpecie operator+(const solidSpecie& st1, const solidSpecie& st2)
{

    scalar sumNmoles = max(st1.nMoles() + st2.nMoles(), SMALL);

    return solidSpecie
    (
        sumNmoles,
        st1.nMoles()/sumNmoles*st1.W()
      + st2.nMoles()/sumNmoles*st2.W()
    );
}


inline solidSpecie operator-(const solidSpecie& st1, const solidSpecie& st2)
{
    scalar diffNmoles = st1.nMoles() - st2.nMoles();
    if (mag(diffNmoles) < SMALL)
    {
        diffNmoles = SMALL;
    }

    return solidSpecie
    (
        diffNmoles,
        st1.nMoles()/diffNmoles*st1.W()
      - st2.nMoles()/diffNmoles*st2.W()
    );
}


inline solidSpecie operator*(const scalar s, const solidSpecie& st)
{
    return solidSpecie
    (
        s*st.nMoles(),
        st.W()
    );
}


inline solidSpecie operator==(const solidSpecie& st1, const solidSpecie& st2)
{
    return st2 - st1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
