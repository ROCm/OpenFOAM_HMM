/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "ensightPart.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightPart, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// TODO - move elsewhere
#if 0
bool Foam::ensightPart::isFieldDefined
(
    const List<scalar>& field
    // const labelUList& addr = cellIds() or faceIds()
) const
{
    forAll(addr, elemI)
    {
        const label id = addr[i];

        if (id >= field.size() || std::isnan(field[id]))
        {
            return false;
        }
    }
    return true;
}
#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightPart::ensightPart(const string& description)
:
    name_(description)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightPart::~ensightPart()
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::ensightGeoFile& Foam::operator<<
(
    ensightGeoFile& os,
    const ensightPart& part
)
{
    part.write(os);
    return os;
}


// ************************************************************************* //
