/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

#include "areaWrite.H"
#include "areaFields.H"
#include "faMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::areaWrite::writeSurface
(
    surfaceWriter& writer,
    const Field<Type>& values,
    const word& fieldName
)
{
    fileName outputName = writer.write(fieldName, values);

    // Case-local file name with "<case>" to make relocatable

    dictionary propsDict;
    propsDict.add
    (
        "file",
        time_.relativePath(outputName, true)
    );
    setProperty(fieldName, propsDict);
}


template<class GeoField>
void Foam::areaWrite::performAction
(
    surfaceWriter& writer,
    const faMesh& areaMesh,
    const IOobjectList& objects
)
{
    wordList fieldNames;
    if (loadFromFiles_)
    {
        fieldNames = objects.sortedNames<GeoField>(fieldSelection_);
    }
    else
    {
        fieldNames = areaMesh.thisDb().sortedNames<GeoField>(fieldSelection_);
    }

    for (const word& fieldName : fieldNames)
    {
        if (verbose_)
        {
            Info<< "write: " << fieldName << endl;
        }

        if (loadFromFiles_)
        {
            const GeoField fld
            (
                IOobject
                (
                    fieldName,
                    time_.timeName(),
                    areaMesh.thisDb(),
                    IOobject::MUST_READ
                ),
                areaMesh
            );

            writeSurface(writer, fld, fieldName);
        }
        else
        {
            writeSurface
            (
                writer,
                areaMesh.thisDb().lookupObject<GeoField>(fieldName),
                fieldName
            );
        }
    }
}


// ************************************************************************* //
