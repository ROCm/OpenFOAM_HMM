/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2012-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::CodedSource<Type>::read(const dictionary& dict)
{
    codedBase::setCodeContext(coeffs_);

    if (!cellSetOption::read(dict))
    {
        return false;
    }

    coeffs_.readEntry("fields", fieldNames_);
    applied_.setSize(fieldNames_.size(), false);

    dict.readCompat<word>("name", {{"redirectType", 1706}}, name_);

    // Code chunks

    codedBase::append("<codeCorrect>");
    {
        const entry& e =
            coeffs_.lookupEntry("codeCorrect", keyType::LITERAL);

        e.readEntry(codeCorrect_);
        dynamicCodeContext::inplaceExpand(codeCorrect_, coeffs_);

        codedBase::append(codeCorrect_);

        dynamicCodeContext::addLineDirective
        (
            codeCorrect_,
            e.startLineNumber(),
            coeffs_
        );
    }

    codedBase::append("<codeAddSup>");
    {
        const entry& e =
            coeffs_.lookupEntry("codeAddSup", keyType::LITERAL);

        e.readEntry(codeAddSup_);
        dynamicCodeContext::inplaceExpand(codeAddSup_, coeffs_);

        codedBase::append(codeAddSup_);

        dynamicCodeContext::addLineDirective
        (
            codeAddSup_,
            e.startLineNumber(),
            coeffs_
        );
    }

    codedBase::append("<codeConstrain>");
    {
        const entry& e =
            coeffs_.lookupEntryCompat
            (
                "codeConstrain",
                {{ "codeSetValue", 1812 }}, keyType::LITERAL
            );

        e.readEntry(codeConstrain_);
        dynamicCodeContext::inplaceExpand(codeConstrain_, coeffs_);

        codedBase::append(codeConstrain_);

        dynamicCodeContext::addLineDirective
        (
            codeConstrain_,
            e.startLineNumber(),
            coeffs_
        );
    }

    return true;
}


// ************************************************************************* //
