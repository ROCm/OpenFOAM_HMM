/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Original code Copyright (C) 2009-2018 Bernhard Gschaider
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "exprMixedFvPatchField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::exprMixedFvPatchField<Type>::setDebug()
{
    if (expressions::patchExprFieldBase::debug_ && !debug)
    {
        debug = 1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::exprMixedFvPatchField<Type>::exprMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    expressions::patchExprFieldBase(true),  // allowGradient
    driver_(this->patch())
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = scalar(1);
}


template<class Type>
Foam::exprMixedFvPatchField<Type>::exprMixedFvPatchField
(
    const exprMixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    expressions::patchExprFieldBase(ptf),
    driver_(this->patch(), ptf.driver_)
{
    setDebug();
    DebugInFunction << nl;
}


template<class Type>
Foam::exprMixedFvPatchField<Type>::exprMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    expressions::patchExprFieldBase(dict, true),
    driver_(this->patch(), dict)
{
    setDebug();
    DebugInFunction << nl;

    // Basic sanity checks
    if (this->valueExpr_.empty() && this->gradExpr_.empty())
    {
        if (this->valueExpr_.empty())
        {
            FatalIOErrorInFunction(dict)
                << "The valueExpr was not defined!" << nl
                << exit(FatalIOError);
        }
        if (this->gradExpr_.empty())
        {
            FatalIOErrorInFunction(dict)
                << "The gradientExpr was not defined!" << nl
                << exit(FatalIOError);
        }
    }


    driver_.readDict(dict);

    // Similar to fvPatchField constructor, which we have bypassed
    dict.readIfPresent("patchType", this->patchType());

    if (dict.found("refValue"))
    {
        this->refValue() = Field<Type>("refValue", dict, p.size());
    }
    else
    {
        this->refValue() = this->patchInternalField();
    }

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );

        if (!dict.found("refValue"))
        {
            // Ensure refValue has a sensible value for the "update" below
            this->refValue() = Field<Type>("value", dict, p.size());
        }
    }
    else
    {
        fvPatchField<Type>::operator=(this->refValue());

        WarningInFunction
            << "No value defined for "
            << this->internalField().name()
            << " on " << this->patch().name() << " therefore using "
            << "the internal field next to the patch"
            << endl;
    }


    if (dict.found("refGradient"))
    {
        this->refGrad() = Field<Type>("refGradient", dict, p.size());
    }
    else
    {
        this->refGrad() = Zero;
    }

    if (dict.found("valueFraction"))
    {
        this->valueFraction() = Field<scalar>("valueFraction", dict, p.size());
    }
    else
    {
        this->valueFraction() = 1;
    }


    if (this->evalOnConstruct_)
    {
        // For potentialFoam or other solvers that don't evaluate
        this->evaluate();
    }
    else
    {
        // Emulate mixedFvPatchField<Type>::evaluate,
        // but avoid our own updateCoeffs
        if (!this->updated())
        {
            this->mixedFvPatchField<Type>::updateCoeffs();
        }

        Field<Type>::operator=
        (
            this->valueFraction()*this->refValue()
          +
            (1.0 - this->valueFraction())*
            (
                this->patchInternalField()
              + this->refGrad()/this->patch().deltaCoeffs()
            )
        );

        fvPatchField<Type>::evaluate();
    }
}


template<class Type>
Foam::exprMixedFvPatchField<Type>::exprMixedFvPatchField
(
    const exprMixedFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    expressions::patchExprFieldBase(ptf),
    driver_(this->patch(), ptf.driver_)
{
    setDebug();
    DebugInFunction << nl;
}


template<class Type>
Foam::exprMixedFvPatchField<Type>::exprMixedFvPatchField
(
    const exprMixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    expressions::patchExprFieldBase(ptf),
    driver_(this->patch(), ptf.driver_)
{
    setDebug();
    DebugInFunction << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::exprMixedFvPatchField<Type>::updateCoeffs()
{
    if (debug)
    {
        InfoInFunction
            << "Value: " << this->valueExpr_ << nl
            << "Gradient: " << this->gradExpr_ << nl
            << "Fraction: " << this->fracExpr_ << nl
            << "Variables: ";

        driver_.writeVariableStrings(Info) << endl;
    }

    if (this->updated())
    {
        return;
    }

    DebugInFunction << " - updating" << nl;


    // Expression evaluation
    {
        driver_.clearVariables();

        if (this->valueExpr_.empty())
        {
            this->refValue() = Zero;
        }
        else
        {
            this->refValue() = driver_.evaluate<Type>(this->valueExpr_);
        }

        bool evalGrad = !this->gradExpr_.empty();

        if (this->fracExpr_.empty() || this->fracExpr_ == "1")
        {
            evalGrad = false;
            this->valueFraction() = scalar(1);
        }
        else if (this->fracExpr_ == "0")
        {
            this->valueFraction() = Zero;
        }
        else
        {
            this->valueFraction() = driver_.evaluate<scalar>(this->fracExpr_);
        }

        if (evalGrad)
        {
            this->refGrad() = driver_.evaluate<Type>(this->gradExpr_);
        }
        else
        {
            this->refGrad() = Zero;
        }
    }

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::exprMixedFvPatchField<Type>::write(Ostream& os) const
{
    mixedFvPatchField<Type>::write(os);
    expressions::patchExprFieldBase::write(os);

    driver_.writeCommon(os, this->debug_ || debug);
}


// ************************************************************************* //
