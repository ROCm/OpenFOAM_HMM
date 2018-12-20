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


void Foam::functionObjects::residuals::createField(const word& fieldName)
{
    if (!writeFields_)
    {
        return;
    }

    const word residualName("initialResidual:" + fieldName);

    if (!mesh_.foundObject<IOField<scalar>>(residualName))
    {
        IOField<scalar>* fieldPtr =
            new IOField<scalar>
            (
                IOobject
                (
                    residualName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                Field<scalar>(mesh_.nCells(), scalar(0))
            );

        fieldPtr->store();
    }
}


void Foam::functionObjects::residuals::writeField(const word& fieldName) const
{
    const word residualName("initialResidual:" + fieldName);

    const auto* residualPtr = mesh_.findObject<IOField<scalar>>(residualName);

    if (residualPtr)
    {
        volScalarField residual
        (
            IOobject
            (
                residualName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            zeroGradientFvPatchField<scalar>::typeName
        );

        residual.primitiveFieldRef() = *residualPtr;
        residual.correctBoundaryConditions();

        residual.write();
    }
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
    fieldSet_(mesh_),
    writeFields_(false),
    initialised_(false)
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

        writeFields_ = dict.lookupOrDefault("writeFields", false);

        return true;
    }

    return false;
}


bool Foam::functionObjects::residuals::execute()
{
    // Note: delaying initialisation until after first iteration so that
    // we can find wildcard fields
    if (!initialised_)
    {
        writeFileHeader(file());

        if (writeFields_)
        {
            for (const word& fieldName : fieldSet_.selection())
            {
                initialiseField<scalar>(fieldName);
                initialiseField<vector>(fieldName);
                initialiseField<sphericalTensor>(fieldName);
                initialiseField<symmTensor>(fieldName);
                initialiseField<tensor>(fieldName);
            }
        }

        initialised_ = true;
    }

    return true;
}


bool Foam::functionObjects::residuals::write()
{
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

    return true;
}


// ************************************************************************* //
