/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    const Field<Type>& value,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    value_(value)
{}


template<class Type>
Foam::Field<Type> Foam::PatchFunction1Types::ConstantField<Type>::getValue
(
    const word& keyword,
    const dictionary& dict,
    const label len
)
{
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
                fld.setSize(len);
                fld = pTraits<Type>(is);
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "expected keyword 'uniform' or 'constant', found "
                    << firstToken.wordToken()
                    << exit(FatalIOError);
            }
        }
        else
        {
            fld.setSize(len);

            is.putBack(firstToken);
            fld = pTraits<Type>(is);
        }
    }
    return fld;
}


template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const polyPatch& pp,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    value_(getValue(entryName, dict, pp.size()))
{}


template<class Type>
Foam::PatchFunction1Types::ConstantField<Type>::ConstantField
(
    const ConstantField<Type>& cnst
)
:
    PatchFunction1<Type>(cnst),
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
    value_(cnst.value_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PatchFunction1Types::ConstantField<Type>::autoMap
(
    const FieldMapper& mapper
)
{
    value_.autoMap(mapper);
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
    //os  << token::SPACE << value_ << token::END_STATEMENT << nl;
    value_.writeEntry(this->name_, os);
}


// ************************************************************************* //
