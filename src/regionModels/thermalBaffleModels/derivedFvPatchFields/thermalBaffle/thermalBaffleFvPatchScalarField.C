/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd
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

#include "thermalBaffleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "emptyPolyPatch.H"
#include "mappedWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(p, iF),
    owner_(false),
    internal_(true),
    baffle_(),
    dict_(),
    extrudeMeshPtr_()
{}


thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const thermalBaffleFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField
    (
        ptf,
        p,
        iF,
        mapper
    ),
    owner_(ptf.owner_),
    internal_(ptf.internal_),
    baffle_(),
    dict_(ptf.dict_),
    extrudeMeshPtr_()
{}


thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(p, iF, dict),
    owner_(false),
    internal_(true),
    baffle_(),
    dict_(dict),
    extrudeMeshPtr_()
{

    const fvMesh& thisMesh = patch().boundaryMesh().mesh();

    typedef regionModels::thermalBaffleModels::thermalBaffleModel baffle;

    word regionName("none");
    dict_.readIfPresent("region", regionName);

    dict_.readIfPresent("internal", internal_);

    const word baffleName("3DBaffle" + regionName);

    if
    (
        !thisMesh.time().foundObject<fvMesh>(regionName)
        && regionName != "none"
    )
    {
        if (!extrudeMeshPtr_)
        {
            createPatchMesh();
        }

        baffle_.reset(baffle::New(thisMesh, dict).ptr());
        owner_ = true;
        baffle_->rename(baffleName);
    }
}


thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const thermalBaffleFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(ptf, iF),
    owner_(ptf.owner_),
    internal_(ptf.internal_),
    baffle_(),
    dict_(ptf.dict_),
    extrudeMeshPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void thermalBaffleFvPatchScalarField::createPatchMesh()
{
    const fvMesh& thisMesh = patch().boundaryMesh().mesh();

    const word regionName(dict_.get<word>("region"));

    List<polyPatch*> regionPatches(3);
    List<word> patchNames(regionPatches.size());
    List<word> patchTypes(regionPatches.size());
    List<dictionary> dicts(regionPatches.size());

    patchNames[bottomPatchID] = word("bottom");
    patchNames[sidePatchID] = word("side");
    patchNames[topPatchID] = word("top");

    patchTypes[bottomPatchID] = mappedWallPolyPatch::typeName;

    if (internal_)
    {
        patchTypes[topPatchID] = mappedWallPolyPatch::typeName;
    }
    else
    {
        patchTypes[topPatchID] = polyPatch::typeName;
    }

    if (dict_.get<bool>("columnCells"))
    {
        patchTypes[sidePatchID] = emptyPolyPatch::typeName;
    }
    else
    {
        patchTypes[sidePatchID] = polyPatch::typeName;
    }

    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch(), dict_);

    const word coupleGroup(mpp.coupleGroup());

    wordList inGroups(1);
    inGroups[0] = coupleGroup;

    // The bottomPatchID is coupled with this patch
    dicts[bottomPatchID].add("coupleGroup", coupleGroup);
    dicts[bottomPatchID].add("inGroups", inGroups);
    dicts[bottomPatchID].add("sampleMode", mpp.sampleModeNames_[mpp.mode()]);
    dicts[bottomPatchID].add("samplePatch", patch().name());
    dicts[bottomPatchID].add("sampleRegion", thisMesh.name());

    // Internal baffle needs a coupled on the topPatchID
    if (internal_)
    {
        const word coupleGroupSlave =
            coupleGroup.substr(0, coupleGroup.find('_')) + "_slave";

        inGroups[0] = coupleGroupSlave;
        dicts[topPatchID].add("coupleGroup", coupleGroupSlave);
        dicts[topPatchID].add("inGroups", inGroups);
        dicts[topPatchID].add("sampleMode", mpp.sampleModeNames_[mpp.mode()]);
    }


    forAll(regionPatches, patchi)
    {
        dictionary&  patchDict = dicts[patchi];
        patchDict.set("nFaces", 0);
        patchDict.set("startFace", 0);

        regionPatches[patchi] = polyPatch::New
        (
            patchTypes[patchi],
            patchNames[patchi],
            dicts[patchi],
            patchi,
            thisMesh.boundaryMesh()
        ).ptr();
    }

    extrudeMeshPtr_.reset
    (
        new extrudePatchMesh
        (
            thisMesh,
            patch(),
            dict_,
            regionName,
            regionPatches
        )
    );
}


void thermalBaffleFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (owner_)
    {
        baffle_->evolve();
    }

    turbulentTemperatureRadCoupledMixedFvPatchScalarField::updateCoeffs();
}


void thermalBaffleFvPatchScalarField::write(Ostream& os) const
{
    turbulentTemperatureRadCoupledMixedFvPatchScalarField::write(os);

    if (owner_)
    {
        os.writeEntry("extrudeModel", dict_.get<word>("extrudeModel"));

        os.writeEntry("nLayers", dict_.get<label>("nLayers"));

        os.writeEntry("expansionRatio", dict_.get<scalar>("expansionRatio"));

        os.writeEntry("columnCells", dict_.get<Switch>("columnCells"));

        const word extrudeModel(dict_.get<word>("extrudeModel") + "Coeffs");

        dict_.subDict(extrudeModel).writeEntry(extrudeModel, os);

        os.writeEntry("region", dict_.get<word>("region"));

        os.writeEntryIfDifferent<bool>("internal", true, internal_);

        os.writeEntry("active", dict_.get<Switch>("active"));

        dict_.subDict("thermoType").writeEntry("thermoType", os);
        dict_.subDict("mixture").writeEntry("mixture", os);
        dict_.subDict("radiation").writeEntry("radiation", os);
   }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    thermalBaffleFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
