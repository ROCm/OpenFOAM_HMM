/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "filmPyrolysisRadiativeCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::filmModelType&
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
filmModel() const
{
    HashTable<const filmModelType*> models
        = db().time().lookupClass<filmModelType>();

    forAllConstIters(models, iter)
    {
        if (iter()->regionMesh().name() == filmRegionName_)
        {
            return *iter();
        }
    }

    DynamicList<word> modelNames;
    forAllConstIters(models, iter)
    {
        modelNames.append(iter()->regionMesh().name());
    }

    FatalErrorInFunction
        << "Unable to locate film region " << filmRegionName_
        << ".  Available regions include: " << modelNames
        << abort(FatalError);

    return **models.begin();
}


const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
pyrolysisModelType&
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
pyrModel() const
{
    HashTable<const pyrolysisModelType*> models =
        db().time().lookupClass<pyrolysisModelType>();

    forAllConstIters(models, iter)
    {
        if (iter()->regionMesh().name() == pyrolysisRegionName_)
        {
            return *iter();
        }
    }

    DynamicList<word> modelNames;
    forAllConstIters(models, iter)
    {
        modelNames.append(iter()->regionMesh().name());
    }

    FatalErrorInFunction
        << "Unable to locate pyrolysis region " << pyrolysisRegionName_
        << ".  Available regions include: " << modelNames
        << abort(FatalError);

    return **models.begin();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase
    (
        patch(),
        "undefined",
        "undefined",
        "undefined-K",
        "undefined-alpha"
    ),
    filmRegionName_("surfaceFilmProperties"),
    pyrolysisRegionName_("pyrolysisProperties"),
    TnbrName_("undefined-Tnbr"),
    qrName_("undefined-qr"),
    convectiveScaling_(1.0),
    filmDeltaDry_(0.0),
    filmDeltaWet_(0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
(
    const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    qrName_(psf.qrName_),
    convectiveScaling_(psf.convectiveScaling_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_)
{}


filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    filmRegionName_
    (
        dict.getOrDefault<word>("filmRegion", "surfaceFilmProperties")
    ),
    pyrolysisRegionName_
    (
        dict.getOrDefault<word>("pyrolysisRegion", "pyrolysisProperties")
    ),
    TnbrName_(dict.lookup("Tnbr")),
    qrName_(dict.lookup("qr")),
    convectiveScaling_(dict.getOrDefault<scalar>("convectiveScaling", 1)),
    filmDeltaDry_(dict.get<scalar>("filmDeltaDry")),
    filmDeltaWet_(dict.get<scalar>("filmDeltaWet"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::
filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
(
    const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    qrName_(psf.qrName_),
    convectiveScaling_(psf.convectiveScaling_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);
}


void filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField& tiptf =
        refCast
        <
            const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
        >(ptf);

    temperatureCoupledBase::rmap(tiptf, addr);
}


void filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const label patchi = patch().index();
    const label nbrPatchi = mpp.samplePolyPatch().index();
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchi];

    scalarField intFld(patchInternalField());

    const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField&
        nbrField =
        refCast
        <
            const filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
        >
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
        );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);

    scalarField& Tp = *this;

    const scalarField K(this->kappa(*this));
    const scalarField nbrK(nbrField.kappa(*this));

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr(nbrK*nbrPatch.deltaCoeffs());
    mpp.distribute(KDeltaNbr);

    scalarField myKDelta(K*patch().deltaCoeffs());

    scalarList Tfilm(patch().size(), Zero);
    scalarList htcwfilm(patch().size(), Zero);
    scalarList filmDelta(patch().size(), Zero);

    const pyrolysisModelType& pyrolysis = pyrModel();
    const filmModelType& film = filmModel();

    // Obtain Rad heat (qr)
    scalarField qr(patch().size(), Zero);

    label coupledPatchi = -1;
    if (pyrolysisRegionName_ == mesh.name())
    {
        coupledPatchi = patchi;
        if (qrName_ != "none")
        {
            qr = nbrPatch.lookupPatchField<volScalarField, scalar>(qrName_);
            mpp.distribute(qr);
        }
    }
    else if (pyrolysis.primaryMesh().name() == mesh.name())
    {
        coupledPatchi = nbrPatch.index();
        if (qrName_ != "none")
        {
            qr = patch().lookupPatchField<volScalarField, scalar>(qrName_);
        }
    }
    else
    {
        FatalErrorInFunction
            << type() << " condition is intended to be applied to either the "
            << "primary or pyrolysis regions only"
            << exit(FatalError);
    }

    const label filmPatchi = pyrolysis.nbrCoupledPatchID(film, coupledPatchi);

    const scalarField htcw(film.htcw().h()().boundaryField()[filmPatchi]);

    // Obtain htcw
    htcwfilm =
        pyrolysis.mapRegionPatchField
        (
            film,
            coupledPatchi,
            filmPatchi,
            htcw,
            true
        );


    // Obtain Tfilm at the boundary through Ts.
    // NOTE: Tf is not good as at the boundary it will retrieve Tp
    Tfilm = film.Ts().boundaryField()[filmPatchi];
    film.toPrimary(filmPatchi, Tfilm);

    // Obtain delta
    filmDelta =
        pyrolysis.mapRegionPatchField<scalar>
        (
            film,
            "deltaf",
            coupledPatchi,
            true
        );

     // Estimate wetness of the film (1: wet , 0: dry)
     scalarField ratio
     (
        min
        (
            max
            (
                (filmDelta - filmDeltaDry_)/(filmDeltaWet_ - filmDeltaDry_),
                scalar(0)
            ),
            scalar(1)
        )
     );

    scalarField qConv(ratio*htcwfilm*(Tfilm - Tp)*convectiveScaling_);

    scalarField qRad((1.0 - ratio)*qr);

    scalarField alpha(KDeltaNbr - (qRad + qConv)/Tp);

    valueFraction() = alpha/(alpha + (1.0 - ratio)*myKDelta);

    refValue() = ratio*Tfilm + (1.0 - ratio)*(KDeltaNbr*nbrIntFld)/alpha;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Qc = gSum(qConv*patch().magSf());
        scalar qr = gSum(qRad*patch().magSf());
        scalar Qt = gSum((qConv + qRad)*patch().magSf());

        Info<< mesh.name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->internalField().name() << " :" << nl
            << "     convective heat[W] : " << Qc << nl
            << "     radiative heat [W] : " << qr << nl
            << "     total heat     [W] : " << Qt << nl
            << "     wall temperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void filmPyrolysisRadiativeCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>
    (
        "filmRegion",
        "surfaceFilmProperties",
        filmRegionName_
    );
    os.writeEntryIfDifferent<word>
    (
        "pyrolysisRegion",
        "pyrolysisProperties",
        pyrolysisRegionName_
    );
    os.writeEntry("Tnbr", TnbrName_);
    os.writeEntry("qr", qrName_);
    os.writeEntry("convectiveScaling", convectiveScaling_);
    os.writeEntry("filmDeltaDry", filmDeltaDry_);
    os.writeEntry("filmDeltaWet", filmDeltaWet_);
    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    filmPyrolysisRadiativeCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
