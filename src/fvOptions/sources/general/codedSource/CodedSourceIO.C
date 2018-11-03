/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "CodedSource.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::CodedSource<Type>::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readEntry("fields", fieldNames_);
        applied_.setSize(fieldNames_.size(), false);

        dict.readCompat<word>("name", {{"redirectType", 1706}}, name_);

        // Code snippets
        {
            const entry& e =
                coeffs_.lookupEntry("codeCorrect", keyType::LITERAL);

            e.readEntry(codeCorrect_);
            stringOps::inplaceTrim(codeCorrect_);
            stringOps::inplaceExpand(codeCorrect_, coeffs_);
            dynamicCodeContext::addLineDirective
            (
                codeCorrect_,
                e.startLineNumber(),
                coeffs_.name()
            );
        }

        {
            const entry& e =
                coeffs_.lookupEntry("codeAddSup", keyType::LITERAL);

            e.readEntry(codeAddSup_);
            stringOps::inplaceTrim(codeAddSup_);
            stringOps::inplaceExpand(codeAddSup_, coeffs_);
            dynamicCodeContext::addLineDirective
            (
                codeAddSup_,
                e.startLineNumber(),
                coeffs_.name()
            );
        }

        {
            const entry& e =
                coeffs_.lookupEntry("codeSetValue", keyType::LITERAL);

            e.readEntry(codeSetValue_);
            stringOps::inplaceTrim(codeSetValue_);
            stringOps::inplaceExpand(codeSetValue_, coeffs_);
            dynamicCodeContext::addLineDirective
            (
                codeSetValue_,
                e.startLineNumber(),
                coeffs_.name()
            );
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
