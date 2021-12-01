/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009-2018 Bernhard Gschaider
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
#include "dictionaryContent.H"

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
    parent_bctype(p, iF),
    expressions::patchExprFieldBase(),
    dict_(),
    driver_(this->patch())
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = scalar(1);
}


template<class Type>
Foam::exprMixedFvPatchField<Type>::exprMixedFvPatchField
(
    const exprMixedFvPatchField<Type>& rhs,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper),
    expressions::patchExprFieldBase(rhs),
    dict_(rhs.dict_),  // Deep copy
    driver_(this->patch(), rhs.driver_, dict_)
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
    parent_bctype(p, iF),
    expressions::patchExprFieldBase
    (
        dict,
        expressions::patchExprFieldBase::expectedTypes::MIXED_TYPE
    ),
    dict_
    (
        // Copy dictionary without "heavy" data chunks
        dictionaryContent::copyDict
        (
            dict,
            wordList(),  // allow
            wordList     // deny
            ({
                "type",  // redundant
                "value", "refValue", "refGradient", "valueFraction"
            })
        )
    ),
    driver_(this->patch(), dict_)
{
    setDebug();
    DebugInFunction << nl;

    // Require one or both of valueExpr, gradientExpr
    if (this->valueExpr_.empty() && this->gradExpr_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "For " << this->internalField().name() << " on "
            << this->patch().name() << nl
            << "Require either or both: valueExpr and gradientExpr" << nl
            << exit(FatalIOError);
    }

    if (this->fracExpr_.empty())
    {
        // No fractionExpr. Expect only one of valueExpr or gradientExpr
        if (!this->valueExpr_.empty() && !this->gradExpr_.empty())
        {
            IOWarningInFunction(dict)
                << "For " << this->internalField().name() << " on "
                << this->patch().name() << nl
                << "Recommend using fractionExpr when specifying both"
                << " valueExpr and gradientExpr. Assuming a value of 1."
                << nl << endl;
        }
    }
    else if (this->fracExpr_ == "0")
    {
        // Gradient only. Expect gradientExpr
        if (this->gradExpr_.empty())
        {
            IOWarningInFunction(dict)
                << "For " << this->internalField().name() << " on "
                << this->patch().name() << nl
                << "Gradient only, but did not specify gradientExpr."
                << nl << endl;
        }
    }
    else if (this->fracExpr_ == "1")
    {
        // Value only. Expect valueExpr
        if (this->valueExpr_.empty())
        {
            IOWarningInFunction(dict)
                << "For " << this->internalField().name() << " on "
                << this->patch().name() << nl
                << "Value only, but did not specify valueExpr."
                << nl << endl;
        }
    }


    driver_.readDict(dict_);

    // Similar to fvPatchField constructor, which we have bypassed
    dict.readIfPresent("patchType", this->patchType());

    bool needsRefValue = true;
    if (dict.found("refValue"))
    {
        needsRefValue = false;
        this->refValue() = Field<Type>("refValue", dict, p.size());
    }

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );

        if (needsRefValue)
        {
            // Ensure refValue has a sensible value for the "update" below
            this->refValue() = static_cast<const Field<Type>&>(*this);
        }
    }
    else
    {
        if (needsRefValue)
        {
            this->refValue() = this->patchInternalField();
        }

        fvPatchField<Type>::operator=(this->refValue());

        #ifdef FULLDEBUG
        WarningInFunction
            << "No value defined for "
            << this->internalField().name() << " on "
            << this->patch().name() << " - using patch internal field" << endl;
        #endif
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
        this->valueFraction() = scalar(1);
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
            this->parent_bctype::updateCoeffs();
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
    const exprMixedFvPatchField<Type>& rhs
)
:
    parent_bctype(rhs),
    expressions::patchExprFieldBase(rhs),
    dict_(rhs.dict_),  // Deep copy
    driver_(this->patch(), rhs.driver_, dict_)
{
    setDebug();
    DebugInFunction << nl;
}


template<class Type>
Foam::exprMixedFvPatchField<Type>::exprMixedFvPatchField
(
    const exprMixedFvPatchField<Type>& rhs,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(rhs, iF),
    expressions::patchExprFieldBase(rhs),
    dict_(rhs.dict_),  // Deep copy
    driver_(this->patch(), rhs.driver_, dict_)
{
    setDebug();
    DebugInFunction << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::exprMixedFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (debug)
    {
        InfoInFunction
            << "Value: " << this->valueExpr_ << nl
            << "Gradient: " << this->gradExpr_ << nl
            << "Fraction: " << this->fracExpr_ << nl
            << "Variables: ";
        driver_.writeVariableStrings(Info) << nl;
        Info<< "... updating" << endl;
    }


    // Expression evaluation
    {
        bool evalValue = (!this->valueExpr_.empty() && this->valueExpr_ != "0");
        bool evalGrad = (!this->gradExpr_.empty() && this->gradExpr_ != "0");
        bool evalFrac = (!this->fracExpr_.empty());
        scalar fraction = 1;

        // Have one or both of valueExpr, gradientExpr (checked in constructor)

        if (this->valueExpr_.empty())
        {
            // No value expression -> gradient only
            fraction = 0;
            evalValue = false;
            evalFrac = false;
        }
        else if (this->gradExpr_.empty())
        {
            // No gradient expression -> value only
            fraction = 1;
            evalGrad = false;
            evalFrac = false;
        }
        else if (this->fracExpr_.empty())
        {
            // No fractionExpr, but has both valueExpr and gradientExpr
            // -> treat as value only (warning in constructor)
            fraction = 1;
            evalGrad = false;
            evalFrac = false;
        }
        else if (this->fracExpr_ == "0")
        {
            // Gradient only
            fraction = 0;
            evalValue = false;
            evalFrac = false;
        }
        else if (this->fracExpr_ == "1")
        {
            // Value only
            fraction = 1;
            evalGrad = false;
            evalFrac = false;
        }


        driver_.clearVariables();

        if (evalValue)
        {
            this->refValue() = driver_.evaluate<Type>(this->valueExpr_);
        }
        else
        {
            this->refValue() = Zero;
        }

        if (evalGrad)
        {
            this->refGrad() = driver_.evaluate<Type>(this->gradExpr_);
        }
        else
        {
            this->refGrad() = Zero;
        }

        if (evalFrac)
        {
            this->valueFraction() = driver_.evaluate<scalar>(this->fracExpr_);
        }
        else
        {
            this->valueFraction() = fraction;
        }
    }

    this->parent_bctype::updateCoeffs();
}


template<class Type>
void Foam::exprMixedFvPatchField<Type>::write(Ostream& os) const
{
    this->parent_bctype::write(os);
    expressions::patchExprFieldBase::write(os);

    driver_.writeCommon(os, this->debug_ || debug);
}


// ************************************************************************* //
