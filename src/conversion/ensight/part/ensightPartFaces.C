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

#include "ensightPartFaces.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightPartFaces, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::ensightPart::localPoints Foam::ensightPartFaces::calcLocalPoints() const
{
    if (contiguousPoints_)
    {
        localPoints ptList;
        ptList.list = identity(points_.size());
        ptList.nPoints = points_.size();
        return ptList;
    }

    localPoints ptList(points_);
    labelList& usedPoints = ptList.list;
    label nPoints = 0;

    // Add all points from faces
    const labelUList& idList = this->faceIds();

    // Add all points from faces
    forAll(idList, i)
    {
        const label id = idList[i] + start_;
        const face& f = faces_[id];

        forAll(f, fp)
        {
            if (usedPoints[f[fp]] == -1)
            {
                usedPoints[f[fp]] = nPoints++;
            }
        }
    }


    // This is not absolutely necessary, but renumber anyhow
    nPoints = 0;
    forAll(usedPoints, ptI)
    {
        if (usedPoints[ptI] > -1)
        {
            usedPoints[ptI] = nPoints++;
        }
    }

    ptList.nPoints = nPoints;
    return ptList;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightPartFaces::ensightPartFaces
(
    label partIndex,
    const string& description,
    const pointField& points,
    const faceList& faces,
    const bool contiguousPoints
)
:
    ensightFaces(partIndex),
    ensightPart(description),
    start_(0),
    patchIndex_(-1),
    faces_(faces),
    points_(points),
    contiguousPoints_(contiguousPoints)
{
    // Classify the face shapes
    classify(faces);
}


Foam::ensightPartFaces::ensightPartFaces
(
    label partIndex,
    const polyMesh& mesh,
    const polyPatch& patch
)
:
    ensightFaces(partIndex),
    ensightPart(patch.name()),
    start_(patch.start()),
    patchIndex_(patch.index()),
    faces_(mesh.faces()),
    points_(mesh.points()),
    contiguousPoints_(false)
{
    // Classify the face shapes
    classify(patch);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightPartFaces::~ensightPartFaces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightPartFaces::writeConnectivity
(
    ensightGeoFile& os,
    const word& key,
    const faceList& faces,
    const labelUList& idList,
    const labelUList& pointMap
) const
{
    if (idList.empty()) return;

    os.writeKeyword(key);
    os.write(idList.size());
    os.newline();

    // Write (polygon) face sizes
    if (key == "nsided")
    {
        // Write the number of points per face
        forAll(idList, i)
        {
            const label id = idList[i] + start_;
            const face& f = faces[id];

            os.write(f.size());
            os.newline();
        }
    }

    // Write the points describing the face
    forAll(idList, i)
    {
        const label id = idList[i] + start_;
        const face& f = faces[id];

        // Convert global -> local index
        // (note: Ensight indices start with 1)
        forAll(f, fp)
        {
            os.write(pointMap[f[fp]] + 1);
        }
        os.newline();
    }
}


void Foam::ensightPartFaces::writeConnectivity
(
    ensightGeoFile& os,
    const word& key,
    const labelUList& idList,
    const labelUList& pointMap
) const
{
    writeConnectivity
    (
        os,
        key,
        faces_,
        idList,
        pointMap
    );
}


void Foam::ensightPartFaces::write
(
    ensightGeoFile& os,
    const pointField& points
) const
{
    if (size())
    {
        const localPoints ptList = calcLocalPoints();
        const labelUList& pointMap = ptList.list;

        os.beginPart(index(), name());
        os.beginCoordinates(ptList.nPoints);

        for (direction cmpt=0; cmpt < point::nComponents; ++cmpt)
        {
            forAll(pointMap, ptI)
            {
                if (pointMap[ptI] > -1)
                {
                    os.write(points[ptI].component(cmpt));
                    os.newline();
                }
            }
        }

        // write part
        for (label typei=0; typei < ensightFaces::nTypes; ++typei)
        {
            const ensightFaces::elemType what = ensightFaces::elemType(typei);

            writeConnectivity
            (
                os,
                ensightFaces::key(what),
                faceIds(what),
                pointMap
            );
        }
    }
}


void Foam::ensightPartFaces::write(ensightGeoFile& os) const
{
    this->write(os, points_);
}


void Foam::ensightPartFaces::writeSummary(Ostream& os) const
{
    os.beginBlock(type());

    os.writeEntry("id",     index()+1); // Ensight starts with 1
    os.writeEntry("name",   name());
    os.writeEntry("start",  start_);
    os.writeEntry("size",   size());

    os.endBlock() << flush;
}


void Foam::ensightPartFaces::dumpInfo(Ostream& os) const
{
    os.beginBlock(type());

    os.writeEntry("id",     index()+1); // Ensight starts with 1
    os.writeEntry("name",   name());
    os.writeEntry("start",  start_);
    os.writeEntry("size",   size());

    for (label typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const ensightFaces::elemType what = ensightFaces::elemType(typei);
        const labelUList& addr = this->faceIds(what);

        os.writeKeyword(ensightFaces::key(what));

        // DIY flat output
        os << addr.size() << '(';
        forAll(addr, i)
        {
            if (i) os << ' ';
            os << addr[i];
        }
        os << ')' << endEntry;
    }

    os.endBlock() << flush;
}


// ************************************************************************* //
