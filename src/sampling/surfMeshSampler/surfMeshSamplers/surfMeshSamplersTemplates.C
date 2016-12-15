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

#include "surfMeshSamplers.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList
Foam::surfMeshSamplers::acceptType() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    return mesh_.names<VolFieldType>(fieldSelection_);
}

#if 0
template<class Type>
Foam::wordList
Foam::surfMeshSamplers::acceptType
(
    const IOobjectList& objects,
    bool fromFiles
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if (fromFiles_)
    {
        // This should actually be in the caller:
        // IOobjectList objects1 = objects.lookup(fieldSelection_);

        return objects.names(VolFieldType::typeName, fieldSelection_);
    }
    else
    {
        return mesh_.names<VolFieldType>(fieldSelection_);
    }
}
#endif


// ************************************************************************* //
