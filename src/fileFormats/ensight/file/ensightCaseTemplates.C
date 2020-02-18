/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "cloud.H"
#include "ensightPTraits.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::ensightFile>
Foam::ensightCase::newData
(
    const word& name,
    const bool isPointData
) const
{
    autoPtr<ensightFile> output;

    if (Pstream::master())
    {
        const ensight::VarName varName(name);

        output = createDataFile(varName);

        // Description
        output().write
        (
            string
            (
                padded(timeIndex_) / varName
              + " <" + pTraits<Type>::typeName + ">"
            )
        );
        output().newline();

        // Remember the field variable for later use
        noteVariable(varName, ensightPTraits<Type>::typeName);

        // Could warn about existing variables that changed representation
        if (isPointData)
        {
            nodeVariables_.set(varName);
        }
    }

    return output;
}


template<class Type>
Foam::autoPtr<Foam::ensightFile>
Foam::ensightCase::newPointData
(
    const word& name
) const
{
    return newData<Type>(name, true);  // POINT_DATA
}


template<class Type>
Foam::autoPtr<Foam::ensightFile>
Foam::ensightCase::newCloudData
(
    const word& cloudName,
    const word& name
) const
{
    autoPtr<Foam::ensightFile> output;

    if (Pstream::master())
    {
        const ensight::VarName varName(name);
        output = createCloudFile(cloudName, varName);

        // Description
        output().write
        (
            string
            (
                padded(timeIndex_) / cloudName / varName
              + " <" + pTraits<Type>::typeName + ">"
            )
        );
        output().newline();

        // Remember the cloud variable for later use
        noteCloud(cloudName, varName, ensightPTraits<Type>::typeName);
    }

    return output;
}


// ************************************************************************* //
