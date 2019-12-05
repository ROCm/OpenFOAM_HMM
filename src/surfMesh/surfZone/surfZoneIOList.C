/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "surfZoneIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfZoneIOList, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfZoneIOList::surfZoneIOList
(
    const IOobject& io
)
:
    regIOobject(io),
    surfZoneList()
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        surfZoneList& zones = *this;

        Istream& is = readStream(typeName);

        PtrList<entry> dictEntries(is);
        zones.setSize(dictEntries.size());

        label facei = 0;
        forAll(zones, zoneI)
        {
            const dictionary& dict = dictEntries[zoneI].dict();

            const label zoneSize = dict.get<label>("nFaces");
            const label startFacei = dict.get<label>("startFace");

            zones[zoneI] = surfZone
            (
                dictEntries[zoneI].keyword(),
                zoneSize,
                startFacei,
                zoneI
            );

            word geoType;
            if (dict.readIfPresent("geometricType", geoType))
            {
                zones[zoneI].geometricType() = geoType;
            }

            if (startFacei != facei)
            {
                FatalErrorInFunction
                    << "surfZones are not ordered. Start of zone " << zoneI
                    << " does not correspond to sum of preceding zones." << nl
                    << "while reading " << io.objectPath() << endl
                    << exit(FatalError);
            }

            facei += zoneSize;
        }

        is.check(FUNCTION_NAME);
        close();
    }
}


Foam::surfZoneIOList::surfZoneIOList
(
    const IOobject& io,
    const UList<surfZone>& content
)
:
    regIOobject(io),
    surfZoneList(content)
{}


Foam::surfZoneIOList::surfZoneIOList
(
    const IOobject& io,
    surfZoneList&& content
)
:
    regIOobject(io),
    surfZoneList(std::move(content))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfZoneIOList::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::surfZoneIOList::operator=(const surfZoneIOList& rhs)
{
    surfZoneList::operator=(rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const surfZoneIOList& list)
{
    os  << list.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    for (const surfZone& item : list)
    {
        item.writeDict(os);
    }

    os  << decrIndent << token::END_LIST;

    return os;
}


// ************************************************************************* //
