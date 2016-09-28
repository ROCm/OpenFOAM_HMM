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

InApplication
    foamToEnsight

Description

SourceFiles
    ensightCloud.C
    ensightCloudTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ensightCloud_H
#define ensightCloud_H

#include "ensightFile.H"
#include "fvMesh.H"
#include "Cloud.H"
#include "IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void ensightParticlePositions
(
    const fvMesh& mesh,
    const fileName& dataDir,
    const label timeIndex,
    const word& cloudName,
    const bool dataExists,
    const IOstream::streamFormat format
);


template<class Type>
void ensightCloudField
(
    const IOobject& fieldObject,
    const fileName& dataDir,
    const label timeIndex,
    const word& cloudName,
    const label cloudNo,
    Ostream& ensightCaseFile,
    const bool dataExists,
    const IOstream::streamFormat format
);


template<class Type>
void writeCloudField
(
    const IOField<Type>& field,
    ensightFile& os
);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightCloudTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
