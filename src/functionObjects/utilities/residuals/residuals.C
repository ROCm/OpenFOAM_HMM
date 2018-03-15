/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

void Foam::functionObjects::residuals::writeFileHeader(Ostream& os)
{
    if (!fieldSet_.updateSelection())
    {
        return;
    }

    if (writtenHeader_)
    {
        writeBreak(file());
    }
    else
    {
        writeHeader(os, "Residuals");
    }

    writeCommented(os, "Time");

    for (const word& fieldName : fieldSet_.selection())
    {
        writeFileHeader<scalar>(os, fieldName);
        writeFileHeader<vector>(os, fieldName);
        writeFileHeader<sphericalTensor>(os, fieldName);
        writeFileHeader<symmTensor>(os, fieldName);
        writeFileHeader<tensor>(os, fieldName);
    }

    os << endl;

    writtenHeader_ = true;
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
    fieldSet_(mesh_)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::residuals::~residuals()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::residuals::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        fieldSet_.read(dict);
        return true;
    }

    return false;
}


bool Foam::functionObjects::residuals::execute()
{
    return true;
}


bool Foam::functionObjects::residuals::write()
{
    if (Pstream::master())
    {
        writeFileHeader(file());

        writeTime(file());

        for (const word& fieldName : fieldSet_.selection())
        {
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
