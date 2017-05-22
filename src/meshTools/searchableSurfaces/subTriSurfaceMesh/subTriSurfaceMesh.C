/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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

#include "subTriSurfaceMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(subTriSurfaceMesh, 0);
addToRunTimeSelectionTable(searchableSurface, subTriSurfaceMesh, dict);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::subTriSurfaceMesh::patchNames(const triSurface& s)
{
    const geometricSurfacePatchList& patches = s.patches();

    wordList names(patches.size());
    forAll(patches, patchI)
    {
        names[patchI] = patches[patchI].name();
    }
    return names;
}


Foam::labelList Foam::subTriSurfaceMesh::selectedRegions
(
    const triSurface& s,
    const wordReList& regionNames
)
{
    const wordList names(patchNames(s));

    labelList regions(names.size());

    label compactI = 0;

    forAll(names, regionI)
    {
        const word& name = names[regionI];

        forAll(regionNames, i)
        {
            if (regionNames[i].match(name))
            {
                regions[compactI++] = regionI;
            }
        }
    }
    regions.setSize(compactI);

    return regions;
}


Foam::triSurface Foam::subTriSurfaceMesh::subset
(
    const IOobject& io,
    const dictionary& dict
)
{
    const word subGeomName(dict.lookup("surface"));

    const triSurfaceMesh& s =
        io.db().lookupObject<triSurfaceMesh>(subGeomName);

    const wordReList regionNames(dict.lookup("patches"));

    labelList regionMap(selectedRegions(s, regionNames));

    if (regionMap.size() == 0)
    {
        FatalIOErrorInFunction(dict)
            << "Found no regions in triSurface matching " << regionNames
            << ". Valid regions are " << patchNames(s)
            << exit(FatalIOError);
    }

    labelList reverseRegionMap(s.patches().size(), -1);
    forAll(regionMap, i)
    {
        reverseRegionMap[regionMap[i]] = i;
    }

    boolList isSelected(s.size(), false);
    forAll(s, triI)
    {
        if (reverseRegionMap[s.triSurface::operator[](triI).region()] != -1)
        {
            isSelected[triI] = true;
        }
    }

    labelList pointMap;
    labelList faceMap;
    return s.subsetMesh(isSelected, pointMap, faceMap);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subTriSurfaceMesh::subTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh(io, subset(io, dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::subTriSurfaceMesh::~subTriSurfaceMesh()
{}


// ************************************************************************* //
