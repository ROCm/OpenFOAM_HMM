/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Description
    Miscellaneous collection of functions and template related to Ensight data

SourceFiles
    ensightOutputFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef ensightOutputSerialCloud_H
#define ensightOutputSerialCloud_H

#include "ensightFile.H"
#include "Cloud.H"
#include "polyMesh.H"
#include "IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
namespace ensightSerialCloud
{

//- Write cloud positions
void writePositions
(
    const polyMesh& mesh,
    const word& cloudName,
    autoPtr<ensightFile> output
);


//- Write cloud field
template<class Type>
bool writeCloudField
(
    const IOField<Type>& field,
    ensightFile& os
);


//- Write cloud field
template<class Type>
bool writeCloudField
(
    const IOobject& fieldObject,
    autoPtr<ensightFile> output
);


} // namespace ensightSerialCloud
} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
    #include "ensightOutputSerialCloudTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
