/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd
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

#include "boundaryRadiationPropertiesFvPatchField.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "fvPatchFieldMapper.H"
#include "radiationModel.H"
#include "absorptionEmissionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::radiation::boundaryRadiationPropertiesFvPatchField::methodType,
        3
    >::names[] =
    {
        "solidRadiation",
        "lookup",
        "model"
    };
}

const Foam::NamedEnum
<
    Foam::radiation::boundaryRadiationPropertiesFvPatchField::methodType,
    3
> Foam::radiation::boundaryRadiationPropertiesFvPatchField::methodTypeNames_;


// * * * * * * * * * * * * * * * * Private functions * * * * * * * * * * * * //

Foam::label
Foam::radiation::boundaryRadiationPropertiesFvPatchField::nbrPatchIndex() const
{
    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    return (mpp.samplePolyPatch().index());
}


const Foam::fvMesh&
Foam::radiation::boundaryRadiationPropertiesFvPatchField::nbrRegion() const
{
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

     return (refCast<const fvMesh>(mpp.sampleMesh()));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::radiation::boundaryRadiationPropertiesFvPatchField::
boundaryRadiationPropertiesFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(p, iF),
    method_(LOOKUP),
    dict_(),
    absorptionEmission_(NULL),
    transmissivity_(NULL)
{}


Foam::radiation::boundaryRadiationPropertiesFvPatchField::
boundaryRadiationPropertiesFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    calculatedFvPatchScalarField(p, iF),
    method_(methodTypeNames_.read(dict.lookup("mode"))),
    dict_(dict),
    absorptionEmission_(NULL),
    transmissivity_(NULL)
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            if (!isA<mappedPatchBase>(p.patch()))
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
            const fvMesh& mesh = this->dimensionedInternalField().mesh();

            //if (dict.found("absorptionEmissionModel"))
            {
                absorptionEmission_.reset
                (
                    absorptionEmissionModel::New(dict, mesh).ptr()
                );
            }

           // if (dict.found("transmissivityModel"))
            {
                transmissivity_.reset
                (
                    transmissivityModel::New(dict, mesh).ptr()
                );
            }
        }
        case LOOKUP:
        {
            //Do nothing
        }
        break;
    }

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );

    }
    else
    {
         fvPatchScalarField::operator=(0.0);
    }
}


Foam::radiation::boundaryRadiationPropertiesFvPatchField::
boundaryRadiationPropertiesFvPatchField
(
    const boundaryRadiationPropertiesFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    calculatedFvPatchScalarField(ptf, p, iF, mapper),
    method_(ptf.method_),
    dict_(ptf.dict_),
    absorptionEmission_(NULL),
    transmissivity_(NULL)
{}


Foam::radiation::boundaryRadiationPropertiesFvPatchField::
boundaryRadiationPropertiesFvPatchField
(
    const boundaryRadiationPropertiesFvPatchField& ptf
)
:
    calculatedFvPatchScalarField(ptf),
    method_(ptf.method_),
    dict_(ptf.dict_),
    absorptionEmission_(NULL),
    transmissivity_(NULL)
{}


Foam::radiation::boundaryRadiationPropertiesFvPatchField::
boundaryRadiationPropertiesFvPatchField
(
    const boundaryRadiationPropertiesFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(ptf, iF),
    method_(ptf.method_),
    dict_(ptf.dict_),
    absorptionEmission_(NULL),
    transmissivity_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::radiation::absorptionEmissionModel&
Foam::radiation::boundaryRadiationPropertiesFvPatchField::
absorptionEmission() const
{
    return absorptionEmission_();
}


const Foam::radiation::transmissivityModel&
Foam::radiation::boundaryRadiationPropertiesFvPatchField::
transmissiveModel() const
{
    return transmissivity_();
}


Foam::tmp<Foam::scalarField> Foam::radiation::
boundaryRadiationPropertiesFvPatchField::
emissivity(const label bandI) const
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
                refCast<const mappedPatchBase>(patch().patch());

            mpp.distribute(emissivity);

            const tmp<scalarField> te(new scalarField(emissivity));

            return te;

        }
        break;

        case LOOKUP:
        {
            tmp<scalarField> e
            (
                 new scalarField("emissivity", dict_, patch().size())
            );

            return e;
        }

        case MODEL:
        {
            const label index = patch().index();

            tmp<scalarField> e
            (
                 new scalarField
                 (
                     absorptionEmission_->e(bandI)().boundaryField()[index]
                 )
            );

            return e;
        }

        default:
        {
            FatalErrorInFunction
                << "Please set 'mode' to one of " << methodTypeNames_.toc()
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}


Foam::tmp<Foam::scalarField> Foam::radiation::
boundaryRadiationPropertiesFvPatchField::absorptivity
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
                refCast<const mappedPatchBase>(patch().patch());

            mpp.distribute(absorp);

            const tmp<scalarField> ta(new scalarField(absorp));

            return ta;

        }
        break;

        case MODEL:
        {
            const label index = patch().index();
            tmp<scalarField> a
            (
                 new scalarField
                 (
                     absorptionEmission_->a(bandI)().boundaryField()[index]
                 )
            );
            return a;
        }

        case LOOKUP:
        {
            tmp<scalarField> a
            (
                 new scalarField("absorptivity", dict_, patch().size())
            );

            return a;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << method_ << endl
                << "Please set 'mode' to one of "
                << methodTypeNames_.toc()
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}


Foam::tmp<Foam::scalarField> Foam::radiation::
boundaryRadiationPropertiesFvPatchField::
transmissivity(const label bandI) const
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
                refCast<const mappedPatchBase>(patch().patch());

            mpp.distribute(trans);

            const tmp<scalarField> tt(new scalarField(trans));

            return tt;

        }
        break;

        case MODEL:
        {
            const label index = patch().index();
            tmp<scalarField> tau
            (
                 new scalarField
                 (
                     transmissivity_->tauEff(bandI)().boundaryField()[index]
                 )
            );
            return tau;
        }

        case LOOKUP:
        {
            tmp<scalarField> tau
            (
                new scalarField
                (
                    "transmissivity", dict_, patch().size()
                )
            );
            return tau;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << method_ << endl
                << "Please set 'mode' to one of "
                << methodTypeNames_.toc()
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}



Foam::tmp<Foam::scalarField> Foam::radiation::
boundaryRadiationPropertiesFvPatchField::
reflectivity(const label bandI) const
{
    const tmp<scalarField> tt = transmissivity(bandI);
    const tmp<scalarField> ta = absorptivity(bandI);

    return (1.0 - tt - ta);
}


void Foam::radiation::boundaryRadiationPropertiesFvPatchField::
write(Ostream& os) const
{
    calculatedFvPatchScalarField::write(os);

    os.writeKeyword("mode") << methodTypeNames_[method_]
        << token::END_STATEMENT << nl;

     switch (method_)
     {
        case MODEL:
        {
            word modelType
            (
                word(dict_.lookup("absorptionEmissionModel"))
            );

            os.writeKeyword("absorptionEmissionModel") << modelType
                << token::END_STATEMENT << nl;

            word modelCoeffs(modelType + word("Coeffs"));
            os.writeKeyword(modelCoeffs);

            dict_.subDict(modelCoeffs).write(os);

            modelType = word(dict_.lookup("transmissivityModel"));

            os.writeKeyword("transmissivityModel") << modelType
                << token::END_STATEMENT << nl;

            modelCoeffs = modelType + word("Coeffs");

            os.writeKeyword(modelCoeffs);

            dict_.subDict(modelCoeffs).write(os);

            break;
        }

        case LOOKUP:
        {
            const scalarField emissivity("emissivity", dict_, patch().size());
            emissivity.writeEntry("emissivity", os);

            const scalarField absorptivity
            (
                "absorptivity", dict_, patch().size()
            );
            absorptivity.writeEntry("absorptivity", os);

            const scalarField transmissivity
            (
                "transmissivity", dict_, patch().size()
            );
            transmissivity.writeEntry("transmissivity", os);

            break;
        }

        case SOLIDRADIATION:
        {
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        boundaryRadiationPropertiesFvPatchField
    );
}
}

// ************************************************************************* //
