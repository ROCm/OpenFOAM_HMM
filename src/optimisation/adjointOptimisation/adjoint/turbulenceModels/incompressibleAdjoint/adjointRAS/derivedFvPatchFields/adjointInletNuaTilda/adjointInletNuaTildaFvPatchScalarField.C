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

#include "adjointInletNuaTildaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointInletNuaTildaFvPatchScalarField::adjointInletNuaTildaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    adjointBoundaryCondition(p, iF, word::null)
{}


adjointInletNuaTildaFvPatchScalarField::adjointInletNuaTildaFvPatchScalarField
(
    const adjointInletNuaTildaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    adjointBoundaryCondition(p, iF, ptf.adjointSolverName_)
{}


adjointInletNuaTildaFvPatchScalarField::adjointInletNuaTildaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    adjointBoundaryCondition(p, iF, dict.get<word>("solverName"))
{}


adjointInletNuaTildaFvPatchScalarField::adjointInletNuaTildaFvPatchScalarField
(
    const adjointInletNuaTildaFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    adjointBoundaryCondition(ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void adjointInletNuaTildaFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==(Zero);

    fixedValueFvPatchScalarField::updateCoeffs();
}


tmp<Field<scalar>> adjointInletNuaTildaFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<scalar>>::New(this->size(), pTraits<scalar>::one);
}


tmp<Field<scalar>> adjointInletNuaTildaFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<scalar>>::New(this->size(), pTraits<scalar>::zero);
}


void adjointInletNuaTildaFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
    os.writeEntry("solverName", adjointSolverName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, adjointInletNuaTildaFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
