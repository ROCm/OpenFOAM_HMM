/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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

#include "adjointOutletNuaTildaFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointOutletNuaTildaFluxFvPatchScalarField::
adjointOutletNuaTildaFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    adjointBoundaryCondition(p, iF, word::null)
{}


adjointOutletNuaTildaFluxFvPatchScalarField::
adjointOutletNuaTildaFluxFvPatchScalarField
(
    const adjointOutletNuaTildaFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    adjointBoundaryCondition(p, iF, ptf.adjointSolverName_)
{}


adjointOutletNuaTildaFluxFvPatchScalarField::
adjointOutletNuaTildaFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    adjointBoundaryCondition(p, iF, dict.get<word>("solverName"))
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
}


adjointOutletNuaTildaFluxFvPatchScalarField::
adjointOutletNuaTildaFluxFvPatchScalarField
(
    const adjointOutletNuaTildaFluxFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    adjointBoundaryCondition(tppsf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void adjointOutletNuaTildaFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator == (scalarField(patch().size(), Zero));

    fixedValueFvPatchScalarField::updateCoeffs();
}


tmp<Field<scalar>>
adjointOutletNuaTildaFluxFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<scalar>>::New(this->size(), pTraits<scalar>::zero);
}


tmp<Field<scalar>>
adjointOutletNuaTildaFluxFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<scalar>>::New(this->size(), pTraits<scalar>::zero);
}


tmp<Field<scalar>>
adjointOutletNuaTildaFluxFvPatchScalarField::gradientBoundaryCoeffs() const
{
    return tmp<Field<scalar>>::New(this->size(), pTraits<scalar>::zero);
}


tmp<Field<scalar>>
adjointOutletNuaTildaFluxFvPatchScalarField::gradientInternalCoeffs() const
{
    return tmp<Field<scalar>>::New(this->size(), pTraits<scalar>::zero);
}


void adjointOutletNuaTildaFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
    os.writeEntry("solverName", adjointSolverName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    adjointOutletNuaTildaFluxFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
