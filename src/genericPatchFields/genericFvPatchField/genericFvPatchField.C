/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "genericFvPatchField.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::genericFvPatchField<Type>::genericFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    calculatedFvPatchField<Type>(p, iF)
{
    FatalErrorInFunction
        << "Trying to construct an genericFvPatchField on patch "
        << this->patch().name()
        << " of field " << this->internalField().name()
        << abort(FatalError);
}


template<class Type>
Foam::genericFvPatchField<Type>::genericFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    calculatedFvPatchField<Type>(p, iF, dict),
    actualTypeName_(dict.lookup("type")),
    dict_(dict)
{
    if (!dict.found("value"))
    {
        FatalIOErrorInFunction(dict)
            << "\n    Cannot find 'value' entry"
            << " on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << nl
            << "    which is required to set the"
               " values of the generic patch field." << nl
            << "    (Actual type " << actualTypeName_ << ")" << nl
            << "\n    Please add the 'value' entry to the write function "
               "of the user-defined boundary-condition\n"
            << exit(FatalIOError);
    }

    for (const entry& dEntry : dict_)
    {
        const keyType& key = dEntry.keyword();

        if (key != "type" && key != "value")
        {
            if
            (
                dEntry.isStream()
             && dEntry.stream().size()
            )
            {
                ITstream& is = dEntry.stream();

                // Read first token
                token firstToken(is);

                if
                (
                    firstToken.isWord()
                 && firstToken.wordToken() == "nonuniform"
                )
                {
                    token fieldToken(is);

                    if (!fieldToken.isCompound())
                    {
                        if
                        (
                            fieldToken.isLabel()
                         && fieldToken.labelToken() == 0
                        )
                        {
                            scalarFields_.insert
                            (
                                dEntry.keyword(),
                                autoPtr<scalarField>::New()
                            );
                        }
                        else
                        {
                            FatalIOErrorInFunction(dict)
                                << "\n    token following 'nonuniform' "
                                   "is not a compound"
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<scalar>>::typeName
                    )
                    {
                        auto fPtr = autoPtr<scalarField>::New();

                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<scalar>>>
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorInFunction(dict)
                                << "\n    size of field " << key
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }

                        scalarFields_.insert(key, fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<vector>>::typeName
                    )
                    {
                        auto fPtr = autoPtr<vectorField>::New();

                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<vector>>>
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorInFunction(dict)
                                << "\n    size of field " << key
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }

                        vectorFields_.insert(key, fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<sphericalTensor>>::typeName
                    )
                    {
                        auto fPtr = autoPtr<sphericalTensorField>::New();

                        fPtr->transfer
                        (
                            dynamicCast
                            <
                                token::Compound<List<sphericalTensor>>
                            >
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorInFunction(dict)
                                << "\n    size of field " << key
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }

                        sphericalTensorFields_.insert(key, fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<symmTensor>>::typeName
                    )
                    {
                        auto fPtr = autoPtr<symmTensorField>::New();

                        fPtr->transfer
                        (
                            dynamicCast
                            <
                                token::Compound<List<symmTensor>>
                            >
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorInFunction(dict)
                                << "\n    size of field " << key
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }

                        symmTensorFields_.insert(key, fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<tensor>>::typeName
                    )
                    {
                        auto fPtr = autoPtr<tensorField>::New();

                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<tensor>>>
                            (
                                fieldToken.transferCompoundToken(is)
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorInFunction(dict)
                                << "\n    size of field " << key
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }

                        tensorFields_.insert(key, fPtr);
                    }
                    else
                    {
                        FatalIOErrorInFunction(dict)
                            << "\n    compound " << fieldToken.compoundToken()
                            << " not supported"
                            << "\n    on patch " << this->patch().name()
                            << " of field "
                            << this->internalField().name()
                            << " in file "
                            << this->internalField().objectPath()
                            << exit(FatalIOError);
                    }
                }
                else if
                (
                    firstToken.isWord()
                 && firstToken.wordToken() == "uniform"
                )
                {
                    token fieldToken(is);

                    if (!fieldToken.isPunctuation())
                    {
                        scalarFields_.insert
                        (
                            key,
                            autoPtr<scalarField>::New
                            (
                                this->size(),
                                fieldToken.number()
                            )
                        );
                    }
                    else
                    {
                        // Read as scalarList.
                        is.putBack(fieldToken);

                        scalarList l(is);

                        if (l.size() == vector::nComponents)
                        {
                            vector vs(l[0], l[1], l[2]);

                            vectorFields_.insert
                            (
                                key,
                                autoPtr<vectorField>::New
                                (
                                    this->size(),
                                    vs
                                )
                            );
                        }
                        else if (l.size() == sphericalTensor::nComponents)
                        {
                            sphericalTensor vs(l[0]);

                            sphericalTensorFields_.insert
                            (
                                key,
                                autoPtr<sphericalTensorField>::New
                                (
                                    this->size(),
                                    vs
                                )
                            );
                        }
                        else if (l.size() == symmTensor::nComponents)
                        {
                            symmTensor vs(l[0], l[1], l[2], l[3], l[4], l[5]);

                            symmTensorFields_.insert
                            (
                                key,
                                autoPtr<symmTensorField>::New
                                (
                                    this->size(),
                                    vs
                                )
                            );
                        }
                        else if (l.size() == tensor::nComponents)
                        {
                            tensor vs
                            (
                                l[0], l[1], l[2],
                                l[3], l[4], l[5],
                                l[6], l[7], l[8]
                            );

                            tensorFields_.insert
                            (
                                key,
                                autoPtr<tensorField>::New
                                (
                                    this->size(),
                                    vs
                                )
                            );
                        }
                        else
                        {
                            FatalIOErrorInFunction(dict)
                                << "\n    unrecognised native type " << l
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->internalField().name()
                                << " in file "
                                << this->internalField().objectPath()
                                << exit(FatalIOError);
                        }
                    }
                }
            }
        }
    }
}


template<class Type>
Foam::genericFvPatchField<Type>::genericFvPatchField
(
    const genericFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    calculatedFvPatchField<Type>(ptf, p, iF, mapper),
    actualTypeName_(ptf.actualTypeName_),
    dict_(ptf.dict_)
{
    forAllConstIter
    (
        HashPtrTable<scalarField>,
        ptf.scalarFields_,
        iter
    )
    {
        scalarFields_.insert
        (
            iter.key(),
            autoPtr<scalarField>::New(*iter(), mapper)
        );
    }

    forAllConstIter
    (
        HashPtrTable<vectorField>,
        ptf.vectorFields_,
        iter
    )
    {
        vectorFields_.insert
        (
            iter.key(),
            autoPtr<vectorField>::New(*iter(), mapper)
        );
    }

    forAllConstIter
    (
        HashPtrTable<sphericalTensorField>,
        ptf.sphericalTensorFields_,
        iter
    )
    {
        sphericalTensorFields_.insert
        (
            iter.key(),
            autoPtr<sphericalTensorField>::New(*iter(), mapper)
        );
    }

    forAllConstIter
    (
        HashPtrTable<symmTensorField>,
        ptf.symmTensorFields_,
        iter
    )
    {
        symmTensorFields_.insert
        (
            iter.key(),
            autoPtr<symmTensorField>::New(*iter(), mapper)
        );
    }

    forAllConstIter
    (
        HashPtrTable<tensorField>,
        ptf.tensorFields_,
        iter
    )
    {
        tensorFields_.insert
        (
            iter.key(),
            autoPtr<tensorField>::New(*iter(), mapper)
        );
    }
}


template<class Type>
Foam::genericFvPatchField<Type>::genericFvPatchField
(
    const genericFvPatchField<Type>& ptf
)
:
    calculatedFvPatchField<Type>(ptf),
    actualTypeName_(ptf.actualTypeName_),
    dict_(ptf.dict_),
    scalarFields_(ptf.scalarFields_),
    vectorFields_(ptf.vectorFields_),
    sphericalTensorFields_(ptf.sphericalTensorFields_),
    symmTensorFields_(ptf.symmTensorFields_),
    tensorFields_(ptf.tensorFields_)
{}


template<class Type>
Foam::genericFvPatchField<Type>::genericFvPatchField
(
    const genericFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    calculatedFvPatchField<Type>(ptf, iF),
    actualTypeName_(ptf.actualTypeName_),
    dict_(ptf.dict_),
    scalarFields_(ptf.scalarFields_),
    vectorFields_(ptf.vectorFields_),
    sphericalTensorFields_(ptf.sphericalTensorFields_),
    symmTensorFields_(ptf.symmTensorFields_),
    tensorFields_(ptf.tensorFields_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::genericFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    calculatedFvPatchField<Type>::autoMap(m);

    forAllIter
    (
        HashPtrTable<scalarField>,
        scalarFields_,
        iter
    )
    {
        iter()->autoMap(m);
    }

    forAllIter
    (
        HashPtrTable<vectorField>,
        vectorFields_,
        iter
    )
    {
        iter()->autoMap(m);
    }

    forAllIter
    (
        HashPtrTable<sphericalTensorField>,
        sphericalTensorFields_,
        iter
    )
    {
        iter()->autoMap(m);
    }

    forAllIter
    (
        HashPtrTable<symmTensorField>,
        symmTensorFields_,
        iter
    )
    {
        iter()->autoMap(m);
    }

    forAllIter
    (
        HashPtrTable<tensorField>,
        tensorFields_,
        iter
    )
    {
        iter()->autoMap(m);
    }
}


template<class Type>
void Foam::genericFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    calculatedFvPatchField<Type>::rmap(ptf, addr);

    const genericFvPatchField<Type>& dptf =
        refCast<const genericFvPatchField<Type>>(ptf);

    forAllIter
    (
        HashPtrTable<scalarField>,
        scalarFields_,
        iter
    )
    {
        HashPtrTable<scalarField>::const_iterator dptfIter =
            dptf.scalarFields_.find(iter.key());

        if (dptfIter != dptf.scalarFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIter
    (
        HashPtrTable<vectorField>,
        vectorFields_,
        iter
    )
    {
        HashPtrTable<vectorField>::const_iterator dptfIter =
            dptf.vectorFields_.find(iter.key());

        if (dptfIter != dptf.vectorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIter
    (
        HashPtrTable<sphericalTensorField>,
        sphericalTensorFields_,
        iter
    )
    {
        HashPtrTable<sphericalTensorField>::const_iterator dptfIter =
            dptf.sphericalTensorFields_.find(iter.key());

        if (dptfIter != dptf.sphericalTensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIter
    (
        HashPtrTable<symmTensorField>,
        symmTensorFields_,
        iter
    )
    {
        HashPtrTable<symmTensorField>::const_iterator dptfIter =
            dptf.symmTensorFields_.find(iter.key());

        if (dptfIter != dptf.symmTensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIter
    (
        HashPtrTable<tensorField>,
        tensorFields_,
        iter
    )
    {
        HashPtrTable<tensorField>::const_iterator dptfIter =
            dptf.tensorFields_.find(iter.key());

        if (dptfIter != dptf.tensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::genericFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorInFunction
        << "cannot be called for a genericFvPatchField"
           " (actual type " << actualTypeName_ << ")"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "generic boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::genericFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorInFunction
        << "cannot be called for a genericFvPatchField"
           " (actual type " << actualTypeName_ << ")"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "generic boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::genericFvPatchField<Type>::gradientInternalCoeffs() const
{
    FatalErrorInFunction
        << "cannot be called for a genericFvPatchField"
           " (actual type " << actualTypeName_ << ")"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "generic boundary condition."
        << abort(FatalError);

    return *this;
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::genericFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    FatalErrorInFunction
        << "cannot be called for a genericFvPatchField"
           " (actual type " << actualTypeName_ << ")"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "generic boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
const Foam::word& Foam::genericFvPatchField<Type>::actualType() const
{
    return actualTypeName_;
}


template<class Type>
void Foam::genericFvPatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", actualTypeName_);

    for (const entry& dEntry : dict_)
    {
        const keyType& key = dEntry.keyword();

        if (key != "type" && key != "value")
        {
            if
            (
                dEntry.isStream()
             && dEntry.stream().size()
             && dEntry.stream()[0].isWord()
             && dEntry.stream()[0].wordToken() == "nonuniform"
            )
            {
                if (scalarFields_.found(key))
                {
                    scalarFields_.find(key)()
                        ->writeEntry(key, os);
                }
                else if (vectorFields_.found(key))
                {
                    vectorFields_.find(key)()
                        ->writeEntry(key, os);
                }
                else if (sphericalTensorFields_.found(key))
                {
                    sphericalTensorFields_.find(key)()
                        ->writeEntry(key, os);
                }
                else if (symmTensorFields_.found(key))
                {
                    symmTensorFields_.find(key)()
                        ->writeEntry(key, os);
                }
                else if (tensorFields_.found(key))
                {
                    tensorFields_.find(key)()
                        ->writeEntry(key, os);
                }
            }
            else
            {
                dEntry.write(os);
            }
        }
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
