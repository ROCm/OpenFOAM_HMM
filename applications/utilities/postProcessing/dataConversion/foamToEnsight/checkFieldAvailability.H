/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Check field availability for last time.
    Done to avoid mapping 'undefined' when a field only exists as time 0.

Requires
    readFields.H (for the checkData function)

\*---------------------------------------------------------------------------*/

// Initially all possible objects that are available at the final time
List<wordHashSet> availableRegionObjectNames(meshes.size());

forAll(meshes, regioni)
{
    const auto& mesh = meshes[regioni];

    IOobjectList objects(mesh, timeDirs.last().name());

    if (!fieldPatterns.empty())
    {
        objects.filterObjects(fieldPatterns);
    }

    // Remove "*_0" restart fields
    objects.prune_0();

    if (!doPointValues)
    {
        // Prune point fields if disabled
        objects.filterClasses
        (
            [](const word& clsName)
            {
                return fieldTypes::point.found(clsName);
            },
            true // prune
        );
    }

    wordList objectNames(objects.sortedNames());

    // Check availability for all times...
    checkData(mesh, timeDirs, objectNames);

    availableRegionObjectNames[regioni] = objectNames;
}


// ************************************************************************* //
