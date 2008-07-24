/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

InClass
    vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamUpdateInformationFields_H
#define vtkPV3FoamUpdateInformationFields_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class patchType, class meshType>
void Foam::vtkPV3Foam::updateInformationFields
(
    vtkDataArraySelection *arraySelection
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInformationFields" << endl;
    }

    // preserve the currently selected values
    const stringList selectedEntries = getSelectedArrayEntries
    (
        arraySelection
    );
    arraySelection->RemoveAllArrays();

    // Search for list of objects for this time
    IOobjectList objects(dbPtr_(), dbPtr_().timeName());
    // Populate the GUI volume/point field arrays

    //- Add volume fields to GUI
/*
    addFields<GeometricField<label, patchType, meshType> >
    (
        arraySelection,
        objects
    );
*/
    addFields<GeometricField<scalar, patchType, meshType> >
    (
        arraySelection,
        objects
    );
    addFields<GeometricField<vector, patchType, meshType> >
    (
        arraySelection,
        objects
    );
    addFields<GeometricField<sphericalTensor, patchType, meshType> >
    (
        arraySelection,
        objects
    );
    addFields<GeometricField<symmTensor, patchType, meshType> >
    (
        arraySelection,
        objects
    );
    addFields<GeometricField<tensor, patchType, meshType> >
    (
        arraySelection,
        objects
    );

    // restore the currently enabled values
    setSelectedArrayEntries
    (
        arraySelection,
        selectedEntries
    );

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateInformationFields" << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
