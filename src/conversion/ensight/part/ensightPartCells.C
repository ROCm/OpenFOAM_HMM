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

#include "ensightPartCells.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightPartCells, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::ensightPart::localPoints Foam::ensightPartCells::calcLocalPoints() const
{
    localPoints ptList(mesh_.points());
    labelList& usedPoints = ptList.list;
    label nPoints = 0;

    // add all points from cells
    const labelUList& idList = this->cellIds();

    forAll(idList, i)
    {
        const label id = idList[i];
        const labelUList& cFaces = mesh_.cells()[id];

        forAll(cFaces, cFacei)
        {
            const face& f = mesh_.faces()[cFaces[cFacei]];

            forAll(f, fp)
            {
                if (usedPoints[f[fp]] == -1)
                {
                    usedPoints[f[fp]] = nPoints++;
                }
            }
        }
    }

    // this is not absolutely necessary, but renumber anyhow
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

Foam::ensightPartCells::ensightPartCells
(
    label partIndex,
    const polyMesh& mesh
)
:
    ensightCells(partIndex),
    ensightPart("cells"),
    mesh_(mesh)
{
    classify(mesh);
}


Foam::ensightPartCells::ensightPartCells
(
    label partIndex,
    const polyMesh& mesh,
    const labelUList& idList
)
:
    ensightCells(partIndex),
    ensightPart("cells"),
    mesh_(mesh)
{
    classify(mesh, idList);
}


Foam::ensightPartCells::ensightPartCells
(
    label partIndex,
    const polyMesh& mesh,
    const cellZone& cZone
)
:
    ensightCells(partIndex),
    ensightPart(cZone.name()),
    mesh_(mesh)
{
    classify(mesh, cZone);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightPartCells::~ensightPartCells()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightPartCells::writeConnectivity
(
    ensightGeoFile& os,
    const word& key,
    const labelUList& idList,
    const labelUList& pointMap
) const
{
    if (idList.empty()) return;

    os.writeKeyword(key);
    os.write(idList.size());
    os.newline();

    // write polyhedral
    if (key == "nfaced")
    {
        const faceList& meshFaces = mesh_.faces();
        const labelUList& owner = mesh_.faceOwner();

        // write the number of faces per element
        forAll(idList, i)
        {
            const label id = idList[i];
            const labelUList& cFace = mesh_.cells()[id];

            os.write(cFace.size());
            os.newline();
        }

        // write the number of points per element face
        forAll(idList, i)
        {
            const label id = idList[i];
            const labelUList& cFace = mesh_.cells()[id];

            forAll(cFace, facei)
            {
                const face& cf = meshFaces[cFace[facei]];

                os.write(cf.size());
                os.newline();
            }
        }

        // write the points describing each element face
        forAll(idList, i)
        {
            const label id = idList[i];
            const labelUList& cFace = mesh_.cells()[id];

            forAll(cFace, cFacei)
            {
                const label faceId = cFace[cFacei];
                const face& cf = meshFaces[faceId];

                // convert global -> local index
                // (note: Ensight indices start with 1)

                // ensight >= 9 needs consistently oriented nfaced cells
                if (id == owner[faceId])
                {
                    forAll(cf, ptI)
                    {
                        os.write(pointMap[cf[ptI]] + 1);
                    }
                }
                else
                {
                    // as per face::reverseFace(), but without copying

                    os.write(pointMap[cf[0]] + 1);
                    for (label ptI = cf.size()-1; ptI > 0; --ptI)
                    {
                        os.write(pointMap[cf[ptI]] + 1);
                    }
                }

                os.newline();
            }
        }
    }
    else
    {
        // write primitive
        const cellShapeList& shapes = mesh_.cellShapes();

        forAll(idList, i)
        {
            const label id = idList[i];
            const cellShape& cellPoints = shapes[id];

            // convert global -> local index
            // (note: Ensight indices start with 1)
            forAll(cellPoints, ptI)
            {
                os.write(pointMap[cellPoints[ptI]] + 1);
            }
            os.newline();
        }
    }
}


void Foam::ensightPartCells::write
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

        // write each element type
        for (label typei=0; typei < ensightCells::nTypes; ++typei)
        {
            const ensightCells::elemType what = ensightCells::elemType(typei);

            writeConnectivity
            (
                os,
                ensightCells::key(what),
                cellIds(what),
                pointMap
            );
        }
    }
}


void Foam::ensightPartCells::write(ensightGeoFile& os) const
{
    this->write(os, mesh_.points());
}


void Foam::ensightPartCells::writeSummary(Ostream& os) const
{
    os.beginBlock(type());

    os.writeEntry("id",     index()+1); // Ensight starts with 1
    os.writeEntry("name",   name());
    os.writeEntry("size",   size());

    os.endBlock() << flush;
}


void Foam::ensightPartCells::dumpInfo(Ostream& os) const
{
    os.beginBlock(type());

    os.writeEntry("id",     index()+1); // Ensight starts with 1
    os.writeEntry("name",   name());
    os.writeEntry("size",   size());

    for (label typei=0; typei < ensightCells::nTypes; ++typei)
    {
        const ensightCells::elemType what = ensightCells::elemType(typei);
        const labelUList& addr = this->cellIds(what);

        os.writeKeyword(ensightCells::key(what));

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
