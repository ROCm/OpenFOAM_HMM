/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "externalFileSource.H"
#include "fam.H"
#include "faScalarMatrix.H"
#include "zeroGradientFaPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{
    defineTypeNameAndDebug(externalFileSource, 0);
    addToRunTimeSelectionTable(option, externalFileSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::externalFileSource::externalFileSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvPatch& p
)
:
    fa::faceSetOption(sourceName, modelType, dict, p),
    fieldName_(dict.get<word>("fieldName")),
    tableName_(dict.get<word>("tableName")),
    pExt_
    (
        IOobject
        (
            "pExt",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("pExt", dimPressure, Zero),
        zeroGradientFaPatchScalarField::typeName
    ),
    value_
    (
        new PatchFunction1Types::MappedFile<scalar>
        (
            p.patch(),
            "uniformValue",
            dict,
            tableName_,          // field table name
            true                 // face values
        )
    ),
    curTimeIndex_(-1)
{
    fieldNames_.resize(1, fieldName_);

    fa::option::resetApplied();

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fa::externalFileSource::addSup
(
    const areaScalarField& solidMass,
    faMatrix<scalar>& eqn,
    const label fieldi
)
{
    const scalar t = mesh().time().value();

    if (isActive() && t > timeStart() && t < (timeStart() + duration()))
    {
        DebugInfo<< name() << ": applying source to " << eqn.psi().name()<<endl;

        if (curTimeIndex_ != mesh().time().timeIndex())
        {
            pExt_.field() = value_->value(t);
            eqn += pExt_/solidMass;
            curTimeIndex_ = mesh().time().timeIndex();
        }
    }
}


bool Foam::fa::externalFileSource::read(const dictionary& dict)
{
    if (fa::option::read(dict))
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
