/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009-2018 Bernhard Gschaider
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

#include "exprFixedValueFvPatchField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::exprFixedValueFvPatchField<Type>::setDebug()
{
    if (expressions::patchExprFieldBase::debug_ && !debug)
    {
        debug = 1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::exprFixedValueFvPatchField<Type>::exprFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    expressions::patchExprFieldBase(),
    driver_(this->patch())
{}


template<class Type>
Foam::exprFixedValueFvPatchField<Type>::exprFixedValueFvPatchField
(
    const exprFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    expressions::patchExprFieldBase(ptf),
    driver_(this->patch(), ptf.driver_)
{
    setDebug();
    DebugInFunction << nl;
}


template<class Type>
Foam::exprFixedValueFvPatchField<Type>::exprFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedValueFvPatchField<Type>(p, iF),
    expressions::patchExprFieldBase
    (
        dict,
        expressions::patchExprFieldBase::expectedTypes::VALUE_TYPE
    ),
    driver_(this->patch(), dict)
{
    setDebug();
    DebugInFunction << nl;

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
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        (*this) == this->patchInternalField();

        #ifdef FULLDEBUG
        WarningInFunction
            << "No value defined for "
            << this->internalField().name() << " on "
            << this->patch().name() << " - using patch internal field" << endl;
        #endif
    }

    if (this->evalOnConstruct_)
    {
        // For potentialFoam or other solvers that don't evaluate
        this->evaluate();
    }
}


template<class Type>
Foam::exprFixedValueFvPatchField<Type>::exprFixedValueFvPatchField
(
    const exprFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    expressions::patchExprFieldBase(ptf),
    driver_(this->patch(), ptf.driver_)
{
    setDebug();
    DebugInFunction << nl;
}


template<class Type>
Foam::exprFixedValueFvPatchField<Type>::exprFixedValueFvPatchField
(
    const exprFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    expressions::patchExprFieldBase(ptf),
    driver_(this->patch(), ptf.driver_)
{
    setDebug();
    DebugInFunction << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::exprFixedValueFvPatchField<Type>::updateCoeffs()
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
            (*this) == driver_.evaluate<Type>(this->valueExpr_);
        }
        else
        {
            (*this) == Zero;
        }
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::exprFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fixedValueFvPatchField<Type>::write(os);
    expressions::patchExprFieldBase::write(os);

    driver_.writeCommon(os, this->debug_ || debug);
}


// ************************************************************************* //
