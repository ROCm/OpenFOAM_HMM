/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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

InClass
    Foam::readFields

Description

SourceFiles
    readFields.C

\*---------------------------------------------------------------------------*/

#ifndef readFields_H
#define readFields_H

#include "fvMeshSubsetProxy.H"
#include "HashSet.H"
#include "PtrList.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Read the fields, optionally subset, and place on the pointer list
template<class GeoField>
label readFields
(
    const fvMeshSubsetProxy& proxy,
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const wordHashSet& selectedFields,
    PtrList<const GeoField>& fields
);


//- Read the fields, and place on the pointer list
template<class GeoField>
label readFields
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const wordHashSet& selectedFields,
    PtrList<const GeoField>& fields
);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "readFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
