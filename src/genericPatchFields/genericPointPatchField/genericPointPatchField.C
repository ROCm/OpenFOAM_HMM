/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "genericPointPatchField.H"
#include "pointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::genericPointPatchField<Type>::genericPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    calculatedPointPatchField<Type>(p, iF)
{
    NotImplemented;
}


template<class Type>
Foam::genericPointPatchField<Type>::genericPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    calculatedPointPatchField<Type>(p, iF, dict),
    actualTypeName_(dict.get<word>("type")),
    dict_(dict)
{
    for (const entry& dEntry : dict_)
    {
        const keyType& key = dEntry.keyword();

        if (key != "type")
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
                                key,
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
            }
        }
    }
}


template<class Type>
Foam::genericPointPatchField<Type>::genericPointPatchField
(
    const genericPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    calculatedPointPatchField<Type>(ptf, p, iF, mapper),
    actualTypeName_(ptf.actualTypeName_),
    dict_(ptf.dict_)
{
    forAllConstIters(ptf.scalarFields_, iter)
    {
        scalarFields_.insert
        (
            iter.key(),
            autoPtr<scalarField>::New(*iter(), mapper)
        );
    }

    forAllConstIters(ptf.vectorFields_, iter)
    {
        vectorFields_.insert
        (
            iter.key(),
            autoPtr<vectorField>::New(*iter(), mapper)
        );
    }

    forAllConstIters(ptf.sphericalTensorFields_, iter)
    {
        sphericalTensorFields_.insert
        (
            iter.key(),
            autoPtr<sphericalTensorField>::New(*iter(), mapper)
        );
    }

    forAllConstIters(ptf.symmTensorFields_, iter)
    {
        symmTensorFields_.insert
        (
            iter.key(),
            autoPtr<symmTensorField>::New(*iter(), mapper)
        );
    }

    forAllConstIters(ptf.tensorFields_, iter)
    {
        tensorFields_.insert
        (
            iter.key(),
            autoPtr<tensorField>::New(*iter(), mapper)
        );
    }
}


template<class Type>
Foam::genericPointPatchField<Type>::genericPointPatchField
(
    const genericPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    calculatedPointPatchField<Type>(ptf, iF),
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
void Foam::genericPointPatchField<Type>::autoMap
(
    const pointPatchFieldMapper& m
)
{
    forAllIters(scalarFields_, iter)
    {
        iter()->autoMap(m);
    }

    forAllIters(vectorFields_, iter)
    {
        iter()->autoMap(m);
    }

    forAllIters(sphericalTensorFields_, iter)
    {
        iter()->autoMap(m);
    }

    forAllIters(symmTensorFields_, iter)
    {
        iter()->autoMap(m);
    }

    forAllIters(tensorFields_, iter)
    {
        iter()->autoMap(m);
    }
}


template<class Type>
void Foam::genericPointPatchField<Type>::rmap
(
    const pointPatchField<Type>& ptf,
    const labelList& addr
)
{
    const genericPointPatchField<Type>& dptf =
        refCast<const genericPointPatchField<Type>>(ptf);

    forAllIters(scalarFields_, iter)
    {
        HashPtrTable<scalarField>::const_iterator dptfIter =
            dptf.scalarFields_.find(iter.key());

        if (dptfIter.found())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIters(vectorFields_, iter)
    {
        HashPtrTable<vectorField>::const_iterator dptfIter =
            dptf.vectorFields_.find(iter.key());

        if (dptfIter.found())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIters(sphericalTensorFields_, iter)
    {
        HashPtrTable<sphericalTensorField>::const_iterator dptfIter =
            dptf.sphericalTensorFields_.find(iter.key());

        if (dptfIter.found())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIters(symmTensorFields_, iter)
    {
        HashPtrTable<symmTensorField>::const_iterator dptfIter =
            dptf.symmTensorFields_.find(iter.key());

        if (dptfIter.found())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    forAllIters(tensorFields_, iter)
    {
        HashPtrTable<tensorField>::const_iterator dptfIter =
            dptf.tensorFields_.find(iter.key());

        if (dptfIter != tensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }
}


template<class Type>
const Foam::word& Foam::genericPointPatchField<Type>::actualType() const
{
    return actualTypeName_;
}


template<class Type>
void Foam::genericPointPatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", actualTypeName_);

    for (const entry& dEntry : dict_)
    {
        const keyType& key = dEntry.keyword();

        if (key != "type")
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
                    scalarFields_.find(key)()->writeEntry(key, os);
                }
                else if (vectorFields_.found(key))
                {
                    vectorFields_.find(key)()->writeEntry(key, os);
                }
                else if (sphericalTensorFields_.found(key))
                {
                    sphericalTensorFields_.find(key)()->writeEntry(key, os);
                }
                else if (symmTensorFields_.found(key))
                {
                    symmTensorFields_.find(key)()->writeEntry(key, os);
                }
                else if (tensorFields_.found(key))
                {
                    tensorFields_.find(key)()->writeEntry(key, os);
                }
            }
            else
            {
                dEntry.write(os);
            }
        }
    }
}


// ************************************************************************* //
