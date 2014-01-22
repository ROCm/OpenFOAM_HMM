/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "BlendedInterfacialModel.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class modelType>
template<class Type>
void Foam::BlendedInterfacialModel<modelType>::correctFixedFluxBCs
(
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    forAll(pair_.phase1().phi().boundaryField(), patchI)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>
            (
                pair_.phase1().phi().boundaryField()[patchI]
            )
        )
        {
            field.boundaryField()[patchI] = pTraits<Type>::zero;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class modelType>
Foam::BlendedInterfacialModel<modelType>::BlendedInterfacialModel
(
    const phasePair::dictTable& modelTable,
    const dictionary& blendingDict,
    const phasePair& pair,
    const orderedPhasePair& pair1In2,
    const orderedPhasePair& pair2In1
)
:
    pair_(pair),
    pair1In2_(pair1In2),
    pair2In1_(pair2In1),
    model_
    (
        modelType::New
        (
            modelTable[pair_],
            pair_
        )
    ),
    model1In2_
    (
        modelType::New
        (
            modelTable[pair1In2_],
            pair1In2_
        )
    ),
    model2In1_
    (
        modelType::New
        (
            modelTable[pair2In1_],
            pair2In1_
        )
    ),
    blending_
    (
        blendingMethod::New
        (
            blendingDict,
            pair1In2_.dispersed(),
            pair1In2_.continuous()
        )
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        blendingDict.lookup("residualAlpha")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class modelType>
Foam::BlendedInterfacialModel<modelType>::~BlendedInterfacialModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class modelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<modelType>::K() const
{
    tmp<volScalarField> f1(blending_->f1());
    tmp<volScalarField> f2(blending_->f2());

    tmp<volScalarField> c
    (
        model_->K()*(f1() - f2())
      + model1In2_->K()*(1 - f1)
      + model2In1_->K()*f2
    );

    c() *= max(pair_.phase1()*pair_.phase2(), residualAlpha_);

    correctFixedFluxBCs(c());

    return c;
}


template<class modelType>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::BlendedInterfacialModel<modelType>::F() const
{
    tmp<volScalarField> f1(blending_->f1());
    tmp<volScalarField> f2(blending_->f2());

    tmp<GeometricField<Type, fvPatchField, volMesh> > v
    (
        model_->F()*(f1() - f2())
      + model1In2_->F()*(1 - f1)
      - model2In1_->F()*f2
    );

    correctFixedFluxBCs(v());

    return v;
}


template<class modelType>
const modelType& Foam::BlendedInterfacialModel<modelType>::phaseModel
(
    const class phaseModel& phase
) const
{
    return &phase == &(pair_.phase1()) ? model1In2_ : model2In1_;
}


// ************************************************************************* //
