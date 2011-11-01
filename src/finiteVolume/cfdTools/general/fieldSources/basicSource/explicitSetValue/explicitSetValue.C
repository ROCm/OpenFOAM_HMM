/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "explicitSetValue.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "HashSet.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(explicitSetValue, 0);
    addToRunTimeSelectionTable
    (
        basicSource,
        explicitSetValue,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::explicitSetValue::setFieldData(const dictionary& dict)
{
    scalarFields_.clear();
    vectorFields_.clear();

    wordList fieldTypes(dict.toc().size());
    wordList fieldNames(dict.toc().size());

    forAll(dict.toc(), i)
    {
        const word& fieldName = dict.toc()[i];
        IOobject io
        (
            fieldName,
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
        if (io.headerOk())
        {
            fieldTypes[i] = io.headerClassName();
            fieldNames[i] = dict.toc()[i];
        }
        else
        {
            FatalErrorIn("explicitSetValue::setFieldData")
                << "header not OK for field " << io.name()
                << abort(FatalError);
        }
    }

    addField(scalarFields_, fieldTypes, fieldNames, dict);
    addField(vectorFields_, fieldTypes, fieldNames, dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::explicitSetValue::explicitSetValue
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    basicSource(name, modelType, dict, mesh)
{
    setFieldData(coeffs_.subDict("fieldData"));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::explicitSetValue::applyToField
(
    const word& fieldName
) const
{
    notImplemented("Foam::explicitSetValue::applyToField(const word&)");
}


void Foam::explicitSetValue::setValue
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    setFieldValue(eqn, scalarFields_[eqn.psi().name()]);
}


void Foam::explicitSetValue::setValue
(
    fvMatrix<vector>& eqn,
    const label
)
{
    setFieldValue(eqn, vectorFields_[eqn.psi().name()]);
}


// ************************************************************************* //
