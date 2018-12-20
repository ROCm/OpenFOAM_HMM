/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::ensightCloud

Description
    A collection of global functions for writing clouds as
    ensight file content.

SourceFiles
    ensightOutputCloud.C
    ensightOutputCloudTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ensightOutputCloud_H
#define ensightOutputCloud_H

#include "ensightFile.H"
#include "autoPtr.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class fvMesh;

namespace ensightCloud
{

//- Write cloud positions
bool writePositions
(
    const fvMesh& mesh,
    const word& cloudName,
    bool exists,
    autoPtr<ensightFile>& output,
    Pstream::commsTypes comm = Pstream::commsTypes::blocking
);


//- Write cloud field, returning true if the field is non-empty.
template<class Type>
bool writeCloudField
(
    const IOField<Type>& field,
    ensightFile& os,
    Pstream::commsTypes comm = Pstream::commsTypes::blocking
);


//- Read cloud field from IOobject (if exists == true) and write,
//- always returning true.
template<class Type>
bool writeCloudField
(
    const IOobject& io,
    bool exists,
    autoPtr<ensightFile>& output,
    Pstream::commsTypes comm = Pstream::commsTypes::blocking
);


} // End namespace ensightCloud
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightOutputCloudTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
