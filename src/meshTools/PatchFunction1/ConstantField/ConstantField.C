/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "ConstantField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const polyPatch& pp,
    const word& entryName,
    const bool isUniform,
    const Type& uniformValue,
    const Field<Type>& nonUniformValue,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    isUniform_(isUniform),
    uniformValue_(uniformValue),
    value_(nonUniformValue)
{
    if (faceValues && nonUniformValue.size() != pp.size())
    {
        FatalIOErrorInFunction(dict)
            << "Supplied field size " << nonUniformValue.size()
            << " is not equal to the number of faces " << pp.size()
            << " of patch " << pp.name() << exit(FatalIOError);
    }
    else if (!faceValues && nonUniformValue.size() != pp.nPoints())
    {
        FatalIOErrorInFunction(dict)
            << "Supplied field size " << nonUniformValue.size()
            << " is not equal to the number of points " << pp.nPoints()
            << " of patch " << pp.name() << exit(FatalIOError);
    }
}


template<class Type>
Foam::Field<Type> Foam::PatchFunction1Types::ConstantField<Type>::getValue
(
    const word& keyword,
    const dictionary& dict,
    const label len,
    bool& isUniform,
    Type& uniformValue
)
{
    isUniform = true;
    uniformValue = Zero;

    Field<Type> fld;

    if (len)
    {
        ITstream& is = dict.lookup(keyword);

        // Read first token
        token firstToken(is);

        if (firstToken.isWord())
        {
            if
            (
                firstToken.wordToken() == "uniform"
             || firstToken.wordToken() == "constant"
            )
            {
                is >> uniformValue;
                fld.setSize(len);
                fld = uniformValue;
            }
            else if (firstToken.wordToken() == "nonuniform")
            {
                List<Type>& list = fld;
                is >> list;
                isUniform = false;

                label currentSize = fld.size();
                if (currentSize != len)
                {
                    if
                    (
                        len < currentSize
                     && FieldBase::allowConstructFromLargerSize
                    )
                    {
                        #ifdef FULLDEBUG
                        IOWarningInFunction(dict)
                            << "Sizes do not match. "
                            << "Re-sizing " << currentSize
                            << " entries to " << len
                            << endl;
                        #endif

                        // Resize the data
                        fld.setSize(len);
                    }
                    else
                    {
                        FatalIOErrorInFunction(dict)
                            << "size " << fld.size()
                            << " is not equal to the given value of " << len
                            << exit(FatalIOError);
                    }
                }
            }
            else
            {
                isUniform = false;
                FatalIOErrorInFunction(dict)
                    << "Expected keyword 'uniform', 'nonuniform' or 'constant'"
                    << ", found " << firstToken.wordToken()
                    << exit(FatalIOError);
            }
        }
        else
        {
            is.putBack(firstToken);
            is >> uniformValue;
            fld.setSize(len);
            fld = uniformValue;
        }
    }

    return fld;
}


template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const polyPatch& pp,
    const word& type,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    value_
    (
        getValue
        (
            entryName,
            dict,
            (faceValues ? pp.size() : pp.nPoints()),
            isUniform_,
            uniformValue_
        )
    )
{}


template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const ConstantField<Type>& cnst
)
:
    PatchFunction1<Type>(cnst),
    isUniform_(cnst.isUniform_),
    uniformValue_(cnst.uniformValue_),
    value_(cnst.value_)
{}


template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const ConstantField<Type>& cnst,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(cnst, pp),
    isUniform_(cnst.isUniform_),
    uniformValue_(cnst.uniformValue_),
    value_(cnst.value_)
{
    // If different sizes do what?
    value_.setSize
    (
        this->faceValues_
      ? this->patch_.size()
      : this->patch_.nPoints()
    );
    if (isUniform_)
    {
        value_ = uniformValue_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PatchFunction1Types::ConstantField<Type>::autoMap
(
    const FieldMapper& mapper
)
{
    value_.autoMap(mapper);

    // If originating from single value override just to make sure
    if (isUniform_)
    {
        value_ = uniformValue_;
    }
}


template<class Type>
void Foam::PatchFunction1Types::ConstantField<Type>::rmap
(
    const PatchFunction1<Type>& pf1,
    const labelList& addr
)
{
    const auto& cst = refCast<const ConstantField<Type>>(pf1);
    value_.rmap(cst.value_, addr);
}


template<class Type>
void Foam::PatchFunction1Types::ConstantField<Type>::writeData
(
    Ostream& os
) const
{
    PatchFunction1<Type>::writeData(os);

    if (isUniform_)
    {
        os.writeKeyword(this->name_)
            << "constant " << uniformValue_
            << token::END_STATEMENT << nl;
    }
    else
    {
        value_.writeEntry(this->name_, os);
    }
}


// ************************************************************************* //
