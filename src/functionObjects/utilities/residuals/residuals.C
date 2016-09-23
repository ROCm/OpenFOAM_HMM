/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "residuals.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(residuals, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        residuals,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::residuals::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Residuals");
    writeCommented(os, "Time");

    forAll(fieldSet_, fieldi)
    {
        const word& fieldName = fieldSet_[fieldi];

        writeFileHeader<scalar>(os, fieldName);
        writeFileHeader<vector>(os, fieldName);
        writeFileHeader<sphericalTensor>(os, fieldName);
        writeFileHeader<symmTensor>(os, fieldName);
        writeFileHeader<tensor>(os, fieldName);
    }

    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::residuals::residuals
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    fieldSet_()
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::residuals::~residuals()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::residuals::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    wordList allFields(dict.lookup("fields"));
    wordHashSet uniqueFields(allFields);
    fieldSet_ = uniqueFields.toc();

    return true;
}


bool Foam::functionObjects::residuals::execute()
{
    return true;
}


bool Foam::functionObjects::residuals::write()
{
    if (Pstream::master())
    {
        writeTime(file());

        forAll(fieldSet_, fieldi)
        {
            const word& fieldName = fieldSet_[fieldi];

            writeResidual<scalar>(fieldName);
            writeResidual<vector>(fieldName);
            writeResidual<sphericalTensor>(fieldName);
            writeResidual<symmTensor>(fieldName);
            writeResidual<tensor>(fieldName);
        }

        file() << endl;
    }

    return true;
}


// ************************************************************************* //
