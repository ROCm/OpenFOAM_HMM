/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::makeSurfaceWriterMethods

Description
    Convenience macros for instantiating writer methods for surfaceWriter
    classes.

\*---------------------------------------------------------------------------*/

#ifndef makeSurfaceWriterMethods_H
#define makeSurfaceWriterMethods_H

#include "surfaceWriter.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSurfaceWriterType(ThisClass)                                       \
    defineTypeNameAndDebug(ThisClass, 0);                                      \
    addToRunTimeSelectionTable(surfaceWriter, ThisClass, word)


#define defineSurfaceWriterWriteField(ThisClass, FieldType)                    \
    Foam::fileName ThisClass::write                                            \
    (                                                                          \
        const fileName& outputDir,                                             \
        const fileName& surfaceName,                                           \
        const meshedSurf& surf,                                                \
        const word& fieldName,                                                 \
        const Field<FieldType>& values,                                        \
        const bool isNodeValues,                                               \
        const bool verbose                                                     \
    ) const                                                                    \
    {                                                                          \
        return writeTemplate                                                   \
        (                                                                      \
            outputDir,                                                         \
            surfaceName,                                                       \
            surf,                                                              \
            fieldName,                                                         \
            values,                                                            \
            isNodeValues,                                                      \
            verbose                                                            \
        );                                                                     \
    }


#define defineSurfaceWriterWriteFields(ThisClass)                              \
    defineSurfaceWriterWriteField(ThisClass, scalar);                          \
    defineSurfaceWriterWriteField(ThisClass, vector);                          \
    defineSurfaceWriterWriteField(ThisClass, sphericalTensor);                 \
    defineSurfaceWriterWriteField(ThisClass, symmTensor);                      \
    defineSurfaceWriterWriteField(ThisClass, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
