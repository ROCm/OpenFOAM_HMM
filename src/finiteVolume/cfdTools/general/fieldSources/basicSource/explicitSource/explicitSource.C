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

#include "explicitSource.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "HashSet.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(explicitSource, 0);
    addToRunTimeSelectionTable
    (
        basicSource,
        explicitSource,
        dictionary
    );


    // * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

    template<> const char* NamedEnum
    <
        explicitSource::volumeModeType,
        2
        >::names[] =
    {
        "absolute",
        "specific"
    };

    const NamedEnum<explicitSource::volumeModeType, 2>
        explicitSource::volumeModeTypeNames_;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::explicitSource::setFieldData(const dictionary& dict)
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
            this->mesh().time().timeName(0),
            this->mesh(),
            IOobject::NO_READ,
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
            FatalErrorIn
            (
                "explicitSource::setFieldData"
            )   << "header not OK " << io.name()
                << exit(FatalError);
        }
    }

    addField(scalarFields_, fieldTypes, fieldNames, dict);
    addField(vectorFields_, fieldTypes, fieldNames, dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::explicitSource::explicitSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    basicSource(name, modelType, dict, mesh),
    dict_(dict.subDict(modelType + "Coeffs")),
    volumeMode_(volumeModeTypeNames_.read(dict_.lookup("volumeMode")))
{
    setFieldData(dict_.subDict("fieldData"));
}


void Foam::explicitSource::addSu(fvMatrix<scalar>& Eqn)
{
    addSource(Eqn, scalarFields_[Eqn.psi().name()]);
}


void Foam::explicitSource::addSu(fvMatrix<vector>& Eqn)
{
    addSource(Eqn, vectorFields_[Eqn.psi().name()]);
}


// ************************************************************************* //
