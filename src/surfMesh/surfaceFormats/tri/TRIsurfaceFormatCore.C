/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "TRIsurfaceFormat.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "StringStream.H"
#include "mergePoints.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::TRIsurfaceFormatCore::TRIsurfaceFormatCore
(
    const fileName& filename
)
:
    sorted_(true),
    points_(0),
    zoneIds_(0),
    sizes_(0)
{
    read(filename);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileFormats::TRIsurfaceFormatCore::~TRIsurfaceFormatCore()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::TRIsurfaceFormatCore::read
(
    const fileName& filename
)
{
    this->clear();
    sorted_ = true;

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    // uses similar structure as STL, just some points
    // the rest of the reader resembles the STL binary reader
    DynamicList<STLpoint> dynPoints;
    DynamicList<label> dynZones;
    DynamicList<label> dynSizes;
    HashTable<label>   lookup;

    // place faces without a group in zone0
    label zoneI = 0;
    dynSizes.append(zoneI);
    lookup.insert("zoneI", zoneI);

    while (is.good())
    {
        string line = this->getLineNoComment(is);

        // handle continuations ?
        //          if (line.back() == '\\')
        //          {
        //              line.resize(line.size()-1);
        //              line += this->getLineNoComment(is);
        //          }

        IStringStream lineStream(line);

        STLpoint p
        (
            readScalar(lineStream),
            readScalar(lineStream),
            readScalar(lineStream)
        );

        if (!lineStream) break;

        dynPoints.append(p);
        dynPoints.append
        (
            STLpoint
            (
                readScalar(lineStream),
                readScalar(lineStream),
                readScalar(lineStream)
            )
        );
        dynPoints.append
        (
            STLpoint
            (
                readScalar(lineStream),
                readScalar(lineStream),
                readScalar(lineStream)
            )
        );

        // zone/colour in .tri file starts with 0x. Skip.
        // ie, instead of having 0xFF, skip 0 and leave xFF to
        // get read as a word and name it "zoneFF"

        char zeroChar;
        lineStream >> zeroChar;

        const word rawName(lineStream);
        const word name("zone" + rawName.substr(1));

        HashTable<label>::const_iterator fnd = lookup.cfind(name);
        if (fnd.found())
        {
            if (zoneI != fnd())
            {
                // group appeared out of order
                sorted_ = false;
            }
            zoneI = fnd();
        }
        else
        {
            zoneI = dynSizes.size();
            lookup.insert(name, zoneI);
            dynSizes.append(0);
        }

        dynZones.append(zoneI);
        dynSizes[zoneI]++;
    }

    // skip empty groups
    label nZone = 0;
    forAll(dynSizes, zoneI)
    {
        if (dynSizes[zoneI])
        {
            if (nZone != zoneI)
            {
                dynSizes[nZone] = dynSizes[zoneI];
            }
            nZone++;
        }
    }
    // truncate addressed size
    dynSizes.setCapacity(nZone);

    // transfer to normal lists
    points_.transfer(dynPoints);
    zoneIds_.transfer(dynZones);
    sizes_.transfer(dynSizes);

    return true;
}


void Foam::fileFormats::TRIsurfaceFormatCore::clear()
{
    sorted_ = true;
    points_.clear();
    zoneIds_.clear();
    sizes_.clear();
}


Foam::label Foam::fileFormats::TRIsurfaceFormatCore::mergePointsMap
(
    labelList& pointMap
) const
{
    // Use merge tolerance as per STL ascii
    return mergePointsMap
    (
        100 * doubleScalarSMALL,
        pointMap
    );
}


Foam::label Foam::fileFormats::TRIsurfaceFormatCore::mergePointsMap
(
    const scalar mergeTol,
    labelList& pointMap
) const
{
    return Foam::mergePoints
    (
        points_,
        mergeTol,
        false, // verbose
        pointMap
    );
}


// ************************************************************************* //
