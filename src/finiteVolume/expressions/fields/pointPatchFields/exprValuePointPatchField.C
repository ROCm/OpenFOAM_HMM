/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
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

#include "exprValuePointPatchField.H"
#include "pointPatchFieldMapper.H"
#include "typeInfo.H"
#include "facePointPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::exprValuePointPatchField<Type>::exprValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    valuePointPatchField<Type>(p, iF),
    expressions::patchExprFieldBase(),
    driver_
    (
        fvPatch::lookupPatch
        (
            dynamicCast<const facePointPatch>(this->patch()).patch()
        )
    )
{}


template<class Type>
Foam::exprValuePointPatchField<Type>::exprValuePointPatchField
(
    const exprValuePointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    valuePointPatchField<Type>(ptf, p, iF, mapper),
    expressions::patchExprFieldBase(ptf),
    driver_
    (
        fvPatch::lookupPatch
        (
            dynamicCast<const facePointPatch>(this->patch()).patch()
        ),
        ptf.driver_
    )
{}


template<class Type>
Foam::exprValuePointPatchField<Type>::exprValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    valuePointPatchField<Type>(p, iF),
    expressions::patchExprFieldBase
    (
        dict,
        expressions::patchExprFieldBase::expectedTypes::VALUE_TYPE,
        true // pointValue
    ),
    driver_
    (
        fvPatch::lookupPatch
        (
            dynamicCast<const facePointPatch>(this->patch()).patch()
        ),
        dict
    )
{
    // Require valueExpr
    if (this->valueExpr_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "The valueExpr was not defined!" << nl
            << exit(FatalIOError);
    }


    driver_.readDict(dict);

    if (dict.found("value"))
    {
        Field<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        WarningInFunction
            << "No value defined for "
            << this->internalField().name()
            << " on " << this->patch().name()
            << endl;
    }

    if (this->evalOnConstruct_)
    {
        // For potentialFoam or other solvers that don't evaluate
        this->evaluate();
    }
}


template<class Type>
Foam::exprValuePointPatchField<Type>::exprValuePointPatchField
(
    const exprValuePointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    valuePointPatchField<Type>(ptf, iF),
    expressions::patchExprFieldBase(ptf),
    driver_
    (
        fvPatch::lookupPatch
        (
            dynamicCast<const facePointPatch>(this->patch()).patch()
        ),
        ptf.driver_
    )
{}


template<class Type>
Foam::exprValuePointPatchField<Type>::exprValuePointPatchField
(
    const exprValuePointPatchField<Type>& ptf
)
:
    valuePointPatchField<Type>(ptf),
    expressions::patchExprFieldBase(ptf),
    driver_
    (
        fvPatch::lookupPatch
        (
            dynamicCast<const facePointPatch>(this->patch()).patch()
        ),
        ptf.driver_
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::exprValuePointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (debug)
    {
        InfoInFunction
            << "Value: " << this->valueExpr_ << nl
            << "Variables: ";
        driver_.writeVariableStrings(Info) << nl;
        Info<< "... updating" << endl;
    }


    // Expression evaluation
    {
        bool evalValue = (!this->valueExpr_.empty() && this->valueExpr_ != "0");


        driver_.clearVariables();

        if (evalValue)
        {
            Field<Type>::operator=
            (
                driver_.evaluate<Type>(this->valueExpr_, true)
            );
        }
        else
        {
            (*this) == Zero;
        }
    }

    valuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::exprValuePointPatchField<Type>::write(Ostream& os) const
{
    valuePointPatchField<Type>::write(os);
    expressions::patchExprFieldBase::write(os);

    this->writeEntry("value", os);

    driver_.writeCommon(os, this->debug_ || debug);
}


// ************************************************************************* //
