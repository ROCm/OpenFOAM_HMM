/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "genericPatchFieldBase.H"
#include "messageStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::genericPatchFieldBase::checkFieldSize
(
    const label fieldSize,
    const label patchSize,
    const word& patchName,
    const keyType& key,
    const IOobject& io
) const
{
    const bool ok = (fieldSize == patchSize);

    if (!ok)
    {
        FatalIOErrorInFunction(dict_)
            << "\n    size of field " << key
            << " (" << fieldSize << ") != patch size (" << patchSize << ')'
            << "\n    on patch " << patchName
            << " of field " << io.name() << " in file "
            << io.objectPath() << nl
            << exit(FatalIOError);
    }

    return ok;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::genericPatchFieldBase::genericPatchFieldBase
(
    const dictionary& dict
)
:
    actualTypeName_(dict.get<word>("type")),
    dict_(dict)
{}


Foam::genericPatchFieldBase::genericPatchFieldBase
(
    const Foam::zero,
    const genericPatchFieldBase& rhs
)
:
    actualTypeName_(rhs.actualTypeName_),
    dict_(rhs.dict_)
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::genericPatchFieldBase::genericFatalSolveError
(
    const word& patchName,
    const IOobject& io
) const
{
    FatalError
        << " (actual type " << actualTypeName_ << ')'
        << "\n    on patch " << patchName
        << " of field " << io.name() << " in file " << io.objectPath() << nl
        << nl
        << "    You are probably trying to solve for a field with a "
           "generic boundary condition." << nl;
}


void Foam::genericPatchFieldBase::reportMissingEntry
(
    const word& entryName,
    const word& patchName,
    const IOobject& io
) const
{
    FatalIOErrorInFunction(dict_)
        << nl
        << "    Missing required '" << entryName << "' entry"
        << " on patch " << patchName
        << " of field " << io.name() << " in file " << io.objectPath() << nl
        << "    (Actual type " << actualTypeName_ << ')' << nl << nl
        << "    Please add the '" << entryName << "' entry to the"
           " write function of the user-defined boundary-condition" << nl
        << exit(FatalIOError);
}


void Foam::genericPatchFieldBase::processGeneric
(
    const label patchSize,
    const word& patchName,
    const IOobject& io,
    const bool separateValue
)
{
    for (const entry& dEntry : dict_)
    {
        const keyType& key = dEntry.keyword();

        if (key == "type" || (separateValue && key == "value"))
        {
            // "type" and possibly "value" handled differently
        }
        else
        {
            processEntry(dEntry, patchSize, patchName, io);
        }
    }
}


bool Foam::genericPatchFieldBase::processEntry
(
    const entry& dEntry,
    const label patchSize,
    const word& patchName,
    const IOobject& io
)
{
    if (!dEntry.isStream())
    {
        return false;
    }

    const keyType& key = dEntry.keyword();
    ITstream& is = dEntry.stream();

    if (is.empty())
    {
        return false;
    }

    #undef FIELDSIZE_CHECK
    #define FIELDSIZE_CHECK(fieldLen)  \
        checkFieldSize(fieldLen, patchSize, patchName, key, io)


    // First token
    token tok(is);

    if (tok.isWord("nonuniform"))
    {
        is >> tok;

        if (tok.isLabel(0))
        {
            // For v2006 and earlier, could have a plain untyped 0
            // without a compound type.
            // Just treat as scalar and hope for the best.
            scalarFields_.insert(key, autoPtr<scalarField>::New());
        }
        else if (!tok.isCompound())
        {
            FatalIOErrorInFunction(dict_)
                << "\n    non-compound token following 'nonuniform'"
                << "\n    on patch " << patchName << " field "
                << io.name() << " in file "
                << io.objectPath() << nl
                << exit(FatalIOError);
            return false;
        }
        else if
        (
            tok.compoundToken().type()
         == token::Compound<List<scalar>>::typeName
        )
        {
            auto fPtr = autoPtr<scalarField>::New();

            fPtr->transfer
            (
                dynamicCast<token::Compound<List<scalar>>>
                (
                    tok.transferCompoundToken(is)
                )
            );

            if (!FIELDSIZE_CHECK(fPtr->size()))
            {
                return false;
            }

            scalarFields_.insert(key, fPtr);
        }
        else if
        (
            tok.compoundToken().type()
         == token::Compound<List<vector>>::typeName
        )
        {
            auto fPtr = autoPtr<vectorField>::New();

            fPtr->transfer
            (
                dynamicCast<token::Compound<List<vector>>>
                (
                    tok.transferCompoundToken(is)
                )
            );

            if (!FIELDSIZE_CHECK(fPtr->size()))
            {
                return false;
            }
            vectorFields_.insert(key, fPtr);
        }
        else if
        (
            tok.compoundToken().type()
         == token::Compound<List<sphericalTensor>>::typeName
        )
        {
            auto fPtr = autoPtr<sphericalTensorField>::New();

            fPtr->transfer
            (
                dynamicCast<token::Compound<List<sphericalTensor>>>
                (
                    tok.transferCompoundToken(is)
                )
            );

            if (!FIELDSIZE_CHECK(fPtr->size()))
            {
                return false;
            }

            sphTensorFields_.insert(key, fPtr);
        }
        else if
        (
            tok.compoundToken().type()
         == token::Compound<List<symmTensor>>::typeName
        )
        {
            auto fPtr = autoPtr<symmTensorField>::New();

            fPtr->transfer
            (
                dynamicCast<token::Compound<List<symmTensor>>>
                (
                    tok.transferCompoundToken(is)
                )
            );

            if (!FIELDSIZE_CHECK(fPtr->size()))
            {
                return false;
            }

            symmTensorFields_.insert(key, fPtr);
        }
        else if
        (
            tok.compoundToken().type()
         == token::Compound<List<tensor>>::typeName
        )
        {
            auto fPtr = autoPtr<tensorField>::New();

            fPtr->transfer
            (
                dynamicCast<token::Compound<List<tensor>>>
                (
                    tok.transferCompoundToken(is)
                )
            );

            if (!FIELDSIZE_CHECK(fPtr->size()))
            {
                return false;
            }

            tensorFields_.insert(key, fPtr);
        }
        else
        {
            FatalIOErrorInFunction(dict_)
                << "\n    unsupported compound " << tok.compoundToken()
                << "\n    on patch " << patchName << " of field "
                << io.name() << " in file "
                << io.objectPath() << nl
                << exit(FatalIOError);
            return false;
        }
    }
    else if (tok.isWord("uniform"))
    {
        is >> tok;

        if (!tok.isPunctuation())
        {
            scalarFields_.insert
            (
                key,
                autoPtr<scalarField>::New(patchSize, tok.number())
            );
        }
        else
        {
            // Read vector-space as list of scalars
            is.putBack(tok);

            scalarList list(is);

            if (list.size() == vector::nComponents)
            {
                vector vs(list[0], list[1], list[2]);

                vectorFields_.insert
                (
                    key,
                    autoPtr<vectorField>::New
                    (
                        patchSize,
                        vs
                    )
                );
            }
            else if (list.size() == sphericalTensor::nComponents)
            {
                sphericalTensor vs(list[0]);

                sphTensorFields_.insert
                (
                    key,
                    autoPtr<sphericalTensorField>::New
                    (
                        patchSize,
                        vs
                    )
                );
            }
            else if (list.size() == symmTensor::nComponents)
            {
                symmTensor vs
                (
                    list[0], list[1], list[2],
                    list[3], list[4],
                    list[5]
                );

                symmTensorFields_.insert
                (
                    key,
                    autoPtr<symmTensorField>::New
                    (
                        patchSize,
                        vs
                    )
                );
            }
            else if (list.size() == tensor::nComponents)
            {
                tensor vs
                (
                    list[0], list[1], list[2],
                    list[3], list[4], list[5],
                    list[6], list[7], list[8]
                );

                tensorFields_.insert
                (
                    key,
                    autoPtr<tensorField>::New
                    (
                        patchSize,
                        vs
                    )
                );
            }
            else
            {
                FatalIOErrorInFunction(dict_)
                    << "\n    unrecognised native type " << flatOutput(list)
                    << "\n    on patch " << patchName << " of field "
                    << io.name() << " in file "
                    << io.objectPath() << nl
                    << exit(FatalIOError);
                return false;
            }
        }
    }

    #undef FIELDSIZE_CHECK

    return true;
}


void Foam::genericPatchFieldBase::putEntry
(
    const entry& e,
    Ostream& os
) const
{
    const keyType& key = e.keyword();

    if
    (
        e.isStream()
     && e.stream().size()
     && e.stream()[0].isWord("nonuniform")
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
        e.write(os);
    }
}


void Foam::genericPatchFieldBase::writeGeneric
(
    Ostream& os,
    const bool separateValue
) const
{
    os.writeEntry("type", actualTypeName_);

    for (const entry& dEntry : dict_)
    {
        const keyType& key = dEntry.keyword();

        if (key == "type" || (separateValue && key == "value"))
        {
            // NB: "type" written first, "value" possibly separately
        }
        else
        {
            putEntry(dEntry, os);
        }
    }
}


void Foam::genericPatchFieldBase::rmapGeneric
(
    const genericPatchFieldBase& rhs,
    const labelList& addr
)
{
    forAllIters(scalarFields_, iter)
    {
        const auto iter2 = rhs.scalarFields_.cfind(iter.key());

        if (iter2.found())
        {
            (*iter)->rmap(*iter2(), addr);
        }
    }

    forAllIters(vectorFields_, iter)
    {
        const auto iter2 = rhs.vectorFields_.cfind(iter.key());

        if (iter2.found())
        {
            (*iter)->rmap(*iter2(), addr);
        }
    }

    forAllIters(sphTensorFields_, iter)
    {
        const auto iter2 = rhs.sphTensorFields_.cfind(iter.key());

        if (iter2.found())
        {
            (*iter)->rmap(*iter2(), addr);
        }
    }

    forAllIters(symmTensorFields_, iter)
    {
        const auto iter2 = rhs.symmTensorFields_.cfind(iter.key());

        if (iter2.found())
        {
            (*iter)->rmap(*iter2(), addr);
        }
    }

    forAllIters(tensorFields_, iter)
    {
        const auto iter2 = rhs.tensorFields_.find(iter.key());

        if (iter2.found())
        {
            (*iter)->rmap(*iter2(), addr);
        }
    }
}


// ************************************************************************* //
