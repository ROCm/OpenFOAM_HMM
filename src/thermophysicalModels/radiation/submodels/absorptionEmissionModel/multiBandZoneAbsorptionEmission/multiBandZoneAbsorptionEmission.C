/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "multiBandZoneAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(multiBandZoneAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            multiBandZoneAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::multiBandZoneAbsorptionEmission::
multiBandZoneAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    absCoeffs_(maxBands_),
    emiCoeffs_(maxBands_),
    nBands_(0),
    zoneAbsorptivity_(),
    zoneEmisivity_(),
    zoneCells_()
{
    coeffsDict_.readEntry("absorptivity", absCoeffs_);
    coeffsDict_.readEntry("emissivity", emiCoeffs_);
    nBands_ = absCoeffs_.size();

    const dictionary& zoneDict = coeffsDict_.subDict("zones");

    zoneDict.readEntry("absorptivity", zoneAbsorptivity_);
    zoneDict.readEntry("emissivity", zoneEmisivity_);

    zoneCells_.setSize(zoneAbsorptivity_.size(), -1);

    label i = 0;
    forAllConstIters(zoneAbsorptivity_, iter)
    {
        label zoneID = mesh.cellZones().findZoneID(iter.key());
        if (zoneID == -1)
        {
            FatalErrorInFunction
                << "Cannot find cellZone " << iter.key() << endl
                << "Valid cellZones are " << mesh.cellZones().names()
                << exit(FatalError);
        }
        zoneCells_[i++] = zoneID;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::multiBandZoneAbsorptionEmission::
~multiBandZoneAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::multiBandZoneAbsorptionEmission::aCont
(
    const label bandI
) const
{
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
            dimensionedScalar("a", dimless/dimLength, absCoeffs_[bandI])
        )
    );

    volScalarField& a = ta.ref();

    forAll(zoneCells_, zonei)
    {
        const cellZone& cZone = mesh().cellZones()[zoneCells_[zonei]];

        tmp<volScalarField> tzoneAbs(a*0.0);
        volScalarField& zoneAbs = tzoneAbs.ref();

        const scalarList& abs = zoneAbsorptivity_.find(cZone.name())();

        forAll(cZone, i)
        {
            label cellId = cZone[i];
            zoneAbs[cellId] = abs[bandI] - absCoeffs_[bandI];
        }

        a += zoneAbs;
    }

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::multiBandZoneAbsorptionEmission::eCont
(
    const label bandI
) const
{
    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("e", dimless/dimLength, emiCoeffs_[bandI])
        )
    );

    volScalarField& e = te.ref();

    forAll(zoneCells_, zonei)
    {
        const cellZone& cZone = mesh().cellZones()[zoneCells_[zonei]];

        tmp<volScalarField> tzoneEm(e*0.0);
        volScalarField& zoneEm = tzoneEm.ref();

        const scalarList& emi = zoneEmisivity_.find(cZone.name())();

        forAll(cZone, i)
        {
            label cellId = cZone[i];
            zoneEm[cellId] = emi[bandI] - emiCoeffs_[bandI];
        }
        e += zoneEm;
    }


    return te;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::multiBandZoneAbsorptionEmission::ECont
(
    const label bandI
) const
{
    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
        )
    );

    return E;
}


// ************************************************************************* //
