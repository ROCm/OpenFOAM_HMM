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

#include "cloud.H"
#include "ensightPTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
Foam::autoPtr<Foam::ensightFile>
Foam::ensightCase::newData
(
    const word& name
) const
{
    autoPtr<ensightFile> output;

    if (Pstream::master())
    {
        const ensight::VarName varName(name);
        output = createDataFile(varName);

        // description
        output().write
        (
            string
            (
                padded(timeIndex_) / varName
              + " <" + pTraits<Type>::typeName + ">"
            )
        );
        output().newline();

        // note variable for later use
        noteVariable(varName, ensightPTraits<Type>::typeName);
    }

    return output;
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

        // description
        output().write
        (
            string
            (
                padded(timeIndex_) / cloudName / varName
              + " <" + pTraits<Type>::typeName + ">"
            )
        );
        output().newline();

        // note cloud variable for later use
        noteCloud(cloudName, varName, ensightPTraits<Type>::typeName);
    }

    return output;
}


// ************************************************************************* //
