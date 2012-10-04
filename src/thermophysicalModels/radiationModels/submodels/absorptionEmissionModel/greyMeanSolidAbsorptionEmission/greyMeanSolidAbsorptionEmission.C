/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "greyMeanSolidAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "zeroGradientFvPatchFields.H"
#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(greyMeanSolidAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            greyMeanSolidAbsorptionEmission,
            dictionary
        );
/*
        template<>
        const char* Foam::NamedEnum
        <
            Foam::radiation::greyMeanSolidAbsorptionEmission::
            radiativePropertiesNumbering,
            3
        >::names[] =
        {
            "absorptivity",
            "emissivity",
            "emission"
        };
*/
    }
}

/*
const Foam::NamedEnum
    <
        Foam::radiation::greyMeanSolidAbsorptionEmission::
        radiativePropertiesNumbering,
        3
    >
    Foam::radiation::greyMeanSolidAbsorptionEmission::
    radiativePropertiesNumberingNames_;
*/

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::radiation::
greyMeanSolidAbsorptionEmission::X
(
    const volScalarField& Yj
) const
{
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    const basicMultiComponentMixture& mixture =
        dynamic_cast<const basicMultiComponentMixture&>(thermo_);

    const label mySpecieI = mixture.species()[Yj.name()];

    tmp<volScalarField> tXj
    (
        new volScalarField
        (
            IOobject
            (
                "Xj",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("Xj", dimless, 0.0),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    scalarField& Xj = tXj().internalField();

    tmp<scalarField> tRhoInv(Xj);
    scalarField& rhoInv = tRhoInv();

    forAll(mixture.Y(), specieI)
    {
        const volScalarField& Yi = mixture.Y(specieI);

        forAll(Xj, iCell)
        {
            rhoInv[iCell] +=
                Yi[iCell]/mixture.rho(specieI, p[iCell], T[iCell]);
        }
    }

    forAll(Xj, iCell)
    {
        Xj[iCell] = Yj[iCell]/mixture.rho(mySpecieI, p[iCell], T[iCell]);
    }

    return (Xj/rhoInv);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyMeanSolidAbsorptionEmission::
greyMeanSolidAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
    thermo_(mesh.lookupObject<fluidThermo>("thermophysicalProperties")),
    speciesNames_(0),
    solidData_
    (
        dynamic_cast
        <
            const basicMultiComponentMixture&
        >(thermo_).species().size()
    )
{
    if (!isA<basicMultiComponentMixture>(thermo_))
    {
        FatalErrorIn
        (
            "radiation::greyMeanSolidAbsorptionEmission::"
            "greyMeanSolidAbsorptionEmission"
            "("
                "const dictionary&, "
                "const fvMesh&"
            ")"
        )   << "Model requires a multi-component thermo package"
            << abort(FatalError);
    }

    label nFunc = 0;
    const dictionary& functionDicts = dict.subDict(typeName + "Coeffs");

    forAllConstIter(dictionary, functionDicts, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        const word& key = iter().keyword();
        speciesNames_.insert(key, nFunc);
        const dictionary& dict = iter().dict();
        dict.lookup("a") >> solidData_[nFunc][absorptivity];
        dict.lookup("e") >> solidData_[nFunc][emissivity];

        nFunc++;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::greyMeanSolidAbsorptionEmission::
~greyMeanSolidAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmission::
calc(const label property) const
{
    const basicMultiComponentMixture& mixture =
        dynamic_cast<const basicMultiComponentMixture&>(thermo_);

    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    scalarField& a = ta().internalField();

    forAllConstIter(HashTable<label>, speciesNames_, iter)
    {
        if (mixture.contains(iter.key()))
        {
            const volScalarField& Y = mixture.Y(iter.key());
            a += solidData_[iter()][property]*X(Y);
        }
    }

    ta().correctBoundaryConditions();
    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmission::eCont
(
    const label bandI
) const
{
   return calc(emissivity);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmission::aCont
(
    const label bandI
) const
{
   return calc(absorptivity);
}

// ************************************************************************* //
