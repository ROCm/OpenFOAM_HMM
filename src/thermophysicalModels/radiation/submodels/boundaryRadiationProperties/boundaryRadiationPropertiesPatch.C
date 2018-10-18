/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "boundaryRadiationPropertiesPatch.H"
#include "mappedPatchBase.H"
#include "radiationModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::radiation::boundaryRadiationPropertiesPatch::methodType
>
Foam::radiation::boundaryRadiationPropertiesPatch::methodTypeNames_
({
    { methodType::SOLIDRADIATION, "solidRadiation" },
    { methodType::LOOKUP, "lookup" },
    { methodType::MODEL, "model" },
});


// * * * * * * * * * * * * * * * * Private functions * * * * * * * * * * * * //

Foam::label
Foam::radiation::boundaryRadiationPropertiesPatch::nbrPatchIndex() const
{
    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch_);

    return (mpp.samplePolyPatch().index());
}


const Foam::fvMesh&
Foam::radiation::boundaryRadiationPropertiesPatch::nbrRegion() const
{
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch_);

     return (refCast<const fvMesh>(mpp.sampleMesh()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::boundaryRadiationPropertiesPatch::
boundaryRadiationPropertiesPatch
(
    const polyPatch& p,
    const dictionary& dict
)
:
    method_(methodTypeNames_.get("mode", dict)),
    dict_(dict),
    absorptionEmission_(nullptr),
    transmissivity_(nullptr),
    patch_(p)
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            if (!isA<mappedPatchBase>(p))
            {
                FatalErrorInFunction
                    << "\n    patch type '" << p.type()
                    << "' not type '" << mappedPatchBase::typeName << "'"
                    << "\n    for patch " << p.name()
                    << abort(FatalIOError);
            }
        }
        break;

        case MODEL:
        {
            const fvMesh& mesh =
                refCast<const fvMesh>(p.boundaryMesh().mesh());

            absorptionEmission_.reset
            (
                absorptionEmissionModel::New(dict, mesh).ptr()
            );

            transmissivity_.reset
            (
                transmissivityModel::New(dict, mesh).ptr()
            );
        }
        case LOOKUP:
        {
            //Do nothing
        }
        break;
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::boundaryRadiationPropertiesPatch::
~boundaryRadiationPropertiesPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::radiation::absorptionEmissionModel&
Foam::radiation::boundaryRadiationPropertiesPatch::absorptionEmission() const
{
    return *absorptionEmission_;
}


const Foam::radiation::transmissivityModel&
Foam::radiation::boundaryRadiationPropertiesPatch::transmissiveModel() const
{
    return *transmissivity_;
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationPropertiesPatch::emissivity
(
    const label bandI
) const
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            const fvMesh& nbrMesh = nbrRegion();

            const radiation::radiationModel& radiation =
                nbrMesh.lookupObject<radiation::radiationModel>
                (
                    "radiationProperties"
                );

            scalarField emissivity
            (
                radiation.absorptionEmission().e(bandI)().boundaryField()
                [
                    nbrPatchIndex()
                ]
            );

            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_);

            mpp.distribute(emissivity);

            return tmp<scalarField>::New(std::move(emissivity));
        }
        break;

        case LOOKUP:
        {
            return tmp<scalarField>::New
            (
                patch_.size(),
                dict_.get<scalar>("emissivity")
            );
        }
        break;

        case MODEL:
        {
            const label index = patch_.index();

            return tmp<scalarField>::New
            (
                absorptionEmission_->e(bandI)().boundaryField()[index]
            );
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "Please set 'mode' to one of "
                << methodTypeNames_
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationPropertiesPatch::absorptivity
(
    const label bandI
) const
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            const fvMesh& nbrMesh = nbrRegion();

            const radiation::radiationModel& radiation =
                nbrMesh.lookupObject<radiation::radiationModel>
                (
                    "radiationProperties"
                );

            scalarField absorp
            (
                radiation.absorptionEmission().a(bandI)().boundaryField()
                [
                    nbrPatchIndex()
                ]
            );

            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_);

            mpp.distribute(absorp);

            return tmp<scalarField>::New(std::move(absorp));

        }
        break;

        case MODEL:
        {
            const label index = patch_.index();

            return tmp<scalarField>::New
            (
                absorptionEmission_->a(bandI)().boundaryField()[index]
            );
        }
        break;

        case LOOKUP:
        {
            return tmp<scalarField>::New
            (
                patch_.size(),
                dict_.get<scalar>("absorptivity")
            );
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << method_ << endl
                << "Please set 'mode' to one of "
                << methodTypeNames_
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationPropertiesPatch::transmissivity
(
    const label bandI
) const
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            const fvMesh& nbrMesh = nbrRegion();

            const radiation::radiationModel& radiation =
                nbrMesh.lookupObject<radiation::radiationModel>
                (
                    "radiationProperties"
                );

            scalarField trans
            (
                radiation.transmissivity().tauEff(bandI)().boundaryField()
                [
                     nbrPatchIndex()
                ]
            );

            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch_);

            mpp.distribute(trans);

            return tmp<scalarField>::New(std::move(trans));

        }
        break;

        case MODEL:
        {
            const label index = patch_.index();

            return tmp<scalarField>::New
            (
                transmissivity_->tauEff(bandI)().boundaryField()[index]
            );
        }

        case LOOKUP:
        {
            return tmp<scalarField>::New
            (
                patch_.size(),
                dict_.get<scalar>("transmissivity")
            );
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << method_ << endl
                << "Please set 'mode' to one of "
                << methodTypeNames_
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationPropertiesPatch::reflectivity
(
    const label bandI
) const
{
    const tmp<scalarField> tt = transmissivity(bandI);
    const tmp<scalarField> ta = absorptivity(bandI);

    return (1.0 - tt - ta);
}


void Foam::radiation::boundaryRadiationPropertiesPatch::write
(
    Ostream& os
) const
{
    os.writeEntry("mode", methodTypeNames_[method_]);

    switch (method_)
    {
        case MODEL:
        {
            word modelType(dict_.get<word>("absorptionEmissionModel"));

            os.writeEntry("absorptionEmissionModel", modelType);

            word modelCoeffs(modelType + word("Coeffs"));
            os.writeKeyword(modelCoeffs);

            dict_.subDict(modelCoeffs).write(os);

            modelType = dict_.get<word>("transmissivityModel");

            os.writeEntry("transmissivityModel", modelType);

            modelCoeffs = modelType + word("Coeffs");

            os.writeKeyword(modelCoeffs);

            dict_.subDict(modelCoeffs).write(os);

            break;
        }

        case LOOKUP:
        {
            const scalarField emissivity("emissivity", dict_, patch_.size());
            emissivity.writeEntry("emissivity", os);

            const scalarField absorptivity
            (
                "absorptivity", dict_, patch_.size()
            );
            absorptivity.writeEntry("absorptivity", os);

            const scalarField transmissivity
            (
                "transmissivity", dict_, patch_.size()
            );
            transmissivity.writeEntry("transmissivity", os);

            break;
        }

        case SOLIDRADIATION:
        {
        }
    }
}


// ************************************************************************* //
