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

#include "surfRegionIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::surfRegionIOList, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfRegionIOList::surfRegionIOList
(
    const IOobject& io
)
:
    surfRegionList(),
    regIOobject(io)
{
    Foam::string functionName =
        "surfRegionIOList::surfRegionIOList"
        "(const IOobject& io)";


    if (readOpt() == IOobject::MUST_READ)
    {
        surfRegionList& regions = *this;

        // read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> dictEntries(is);
        regions.setSize(dictEntries.size());

        label faceI = 0;
        forAll(regions, regionI)
        {
            const dictionary& dict = dictEntries[regionI].dict();

            label regionSize = readLabel(dict.lookup("nFaces"));
            label startFaceI = readLabel(dict.lookup("startFace"));

            regions[regionI] = surfRegion
            (
                dictEntries[regionI].keyword(),
                regionSize,
                startFaceI,
                regionI
            );

            word geoType;
            if (dict.readIfPresent("geometricType", geoType))
            {
                regions[regionI].geometricType() = geoType;
            }

            if (startFaceI != faceI)
            {
                FatalErrorIn(functionName)
                    << "Regions are not ordered. Start of region " << regionI
                    << " does not correspond to sum of preceding regions."
                    << endl
                    << "while reading " << io.objectPath()
                    << exit(FatalError);
            }

            faceI += regionSize;
        }

        // Check state of IOstream
        is.check(functionName.c_str());

        close();
    }
}

// Construct from IOObject
Foam::surfRegionIOList::surfRegionIOList
(
    const IOobject& io,
    const surfRegionList& regions
)
:
    surfRegionList(regions),
    regIOobject(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfRegionIOList::~surfRegionIOList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// writeData member function required by regIOobject
bool Foam::surfRegionIOList::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const surfRegionIOList& L)
{
    os  << L.size() << nl << token::BEGIN_LIST;

    forAll(L, i)
    {
        L[i].writeDict(os);
    }

    os  << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
