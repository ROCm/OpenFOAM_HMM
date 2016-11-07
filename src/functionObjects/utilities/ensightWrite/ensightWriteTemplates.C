/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "Time.H"
#include "ensightOutput.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
int Foam::functionObjects::ensightWrite::writeVolField
(
    const word& inputName,
    int& state
)
{
    // State: return 0 (not-processed), -1 (skip), +1 ok
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    // Already done, or not available
    if (state || !foundObject<VolFieldType>(inputName))
    {
        return state;
    }

    autoPtr<ensightFile> os = ensCase().newData<Type>(inputName);
    ensightOutput::writeField<Type>
    (
        lookupObject<VolFieldType>(inputName),
        ensMesh(),
        os
    );

    Log << " " << inputName;

    state = +1;
    return state;
}


// ************************************************************************* //
