/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
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

Namespace
    Foam::ensightOutput

Description
    A collection of functions for writing ensight file content.

SourceFiles
    ensightOutputTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ensightOutput_H
#define ensightOutput_H

#include "ensightFile.H"
#include "ensightCells.H"
#include "ensightFaces.H"
#include "ensightPTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
namespace ensightOutput
{

/*---------------------------------------------------------------------------*\
                        Namespace ensightOutput::Detail
\*---------------------------------------------------------------------------*/

//- \brief Implementation details and output backends that would not normally
//- be called directly by a user.

namespace Detail
{

//- Write field content (component-wise) for the given ensight element type
template<template<typename> class FieldContainer, class Type>
bool writeFieldComponents
(
    const char* key,
    const FieldContainer<Type>& fld,
    ensightFile& os,
    bool parallel         //!< Collective write?
);


//- Write a field of faces values as an indirect list,
//- using the face ids from ensightFaces
template<class Type>
bool writeFaceField
(
    const Field<Type>& fld,
    const ensightFaces& part,
    ensightFile& os,
    bool parallel         //!< Collective write?
);


//- Write a sub-field of faces values as an indirect list,
//- using the sublist sizing information from ensightFaces
template<class Type>
bool writeFaceSubField
(
    const Field<Type>& fld,
    const ensightFaces& part,
    ensightFile& os,
    bool parallel         //!< Collective write?
);



//- Write a field of cell values as an indirect list,
//- using the cell ids from ensightCells
template<class Type>
bool writeCellField
(
    const Field<Type>& fld,
    const ensightCells& part,
    ensightFile& os,
    bool parallel         //!< Collective write?
);

} // End namespace Detail


/*---------------------------------------------------------------------------*\
                        Namespace ensightOutput::Serial
\*---------------------------------------------------------------------------*/

//- \brief Output routines that are either designed for serial-only,
//- or for which parallelization is pending.
namespace Serial
{

//- Write a field of point (node) values (already compacted?)
template<class Type>
bool writePointField
(
    const Field<Type>& fld,
    const ensightFaces& part,
    ensightFile& os
);

} // End namespace Serial

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ensightOutput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightOutputTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
