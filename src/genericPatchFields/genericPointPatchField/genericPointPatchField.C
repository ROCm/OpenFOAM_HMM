/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    const label patchSize = this->size();

    for (const entry& dEntry : dict_)
    {
        const keyType& key = dEntry.keyword();

        if
        (
            key == "type"
         || !dEntry.isStream() || dEntry.stream().empty()
        )
        {
            continue;
        }


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
                    scalarFields_.insert(key, autoPtr<scalarField>::New());
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
                        << this->internalField().objectPath() << nl
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

                if (fPtr->size() != patchSize)
                {
                    FatalIOErrorInFunction(dict)
                        << "\n    size of field " << key
                        << " (" << fPtr->size() << ')'
                        << " is not the same size as the patch ("
                        << patchSize << ')'
                        << "\n    on patch " << this->patch().name()
                        << " of field "
                        << this->internalField().name()
                        << " in file "
                        << this->internalField().objectPath() << nl
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

                if (fPtr->size() != patchSize)
                {
                    FatalIOErrorInFunction(dict)
                        << "\n    size of field " << key
                        << " (" << fPtr->size() << ')'
                        << " is not the same size as the patch ("
                        << patchSize << ')'
                        << "\n    on patch " << this->patch().name()
                        << " of field "
                        << this->internalField().name()
                        << " in file "
                        << this->internalField().objectPath() << nl
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
                    dynamicCast<token::Compound<List<sphericalTensor>>>
                    (
                        fieldToken.transferCompoundToken(is)
                    )
                );

                if (fPtr->size() != patchSize)
                {
                    FatalIOErrorInFunction(dict)
                        << "\n    size of field " << key
                        << " (" << fPtr->size() << ')'
                        << " is not the same size as the patch ("
                        << patchSize << ')'
                        << "\n    on patch " << this->patch().name()
                        << " of field "
                        << this->internalField().name()
                        << " in file "
                        << this->internalField().objectPath() << nl
                        << exit(FatalIOError);
                }

                sphTensorFields_.insert(key, fPtr);
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
                    dynamicCast<token::Compound<List<symmTensor>>>
                    (
                        fieldToken.transferCompoundToken(is)
                    )
                );

                if (fPtr->size() != patchSize)
                {
                    FatalIOErrorInFunction(dict)
                        << "\n    size of field " << key
                        << " (" << fPtr->size() << ')'
                        << " is not the same size as the patch ("
                        << patchSize << ')'
                        << "\n    on patch " << this->patch().name()
                        << " of field "
                        << this->internalField().name()
                        << " in file "
                        << this->internalField().objectPath() << nl
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

                if (fPtr->size() != patchSize)
                {
                    FatalIOErrorInFunction(dict)
                        << "\n    size of field " << key
                        << " (" << fPtr->size() << ')'
                        << " is not the same size as the patch ("
                        << patchSize << ')'
                        << "\n    on patch " << this->patch().name()
                        << " of field "
                        << this->internalField().name()
                        << " in file "
                        << this->internalField().objectPath() << nl
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
                    << this->internalField().objectPath() << nl
                    << exit(FatalIOError);
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

    forAllConstIters(ptf.sphTensorFields_, iter)
    {
        sphTensorFields_.insert
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
    sphTensorFields_(ptf.sphTensorFields_),
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
        (*iter)->autoMap(m);
    }

    forAllIters(vectorFields_, iter)
    {
        (*iter)->autoMap(m);
    }

    forAllIters(sphTensorFields_, iter)
    {
        (*iter)->autoMap(m);
    }

    forAllIters(symmTensorFields_, iter)
    {
        (*iter)->autoMap(m);
    }

    forAllIters(tensorFields_, iter)
    {
        (*iter)->autoMap(m);
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
        const auto iter2 = dptf.scalarFields_.cfind(iter.key());

        if (iter2.found())
        {
            (*iter)->rmap(*iter2(), addr);
        }
    }

    forAllIters(vectorFields_, iter)
    {
        const auto iter2 = dptf.vectorFields_.cfind(iter.key());

        if (iter2.found())
        {
            (*iter)->rmap(*iter2(), addr);
        }
    }

    forAllIters(sphTensorFields_, iter)
    {
        const auto iter2 = dptf.sphTensorFields_.find(iter.key());

        if (iter2.found())
        {
            (*iter)->rmap(*iter2(), addr);
        }
    }

    forAllIters(symmTensorFields_, iter)
    {
        const auto iter2 = dptf.symmTensorFields_.find(iter.key());

        if (iter2.found())
        {
            (*iter)->rmap(*iter2(), addr);
        }
    }

    forAllIters(tensorFields_, iter)
    {
        const auto iter2 = dptf.tensorFields_.find(iter.key());

        if (iter2.found())
        {
            (*iter)->rmap(*iter2(), addr);
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

        if (key == "type")
        {
            // NB: "type" written first, no special treatment for "value"
            continue;
        }
        else if
        (
            dEntry.isStream()
         && dEntry.stream().size()
         && dEntry.stream()[0].isWord()
         && dEntry.stream()[0].wordToken() == "nonuniform"
        )
        {
            if (scalarFields_.found(key))
            {
                scalarFields_.cfind(key)()->writeEntry(key, os);
            }
            else if (vectorFields_.found(key))
            {
                vectorFields_.cfind(key)()->writeEntry(key, os);
            }
            else if (sphTensorFields_.found(key))
            {
                sphTensorFields_.cfind(key)()->writeEntry(key, os);
            }
            else if (symmTensorFields_.found(key))
            {
                symmTensorFields_.cfind(key)()->writeEntry(key, os);
            }
            else if (tensorFields_.found(key))
            {
                tensorFields_.cfind(key)()->writeEntry(key, os);
            }
        }
        else
        {
            dEntry.write(os);
        }
    }
}


// ************************************************************************* //
