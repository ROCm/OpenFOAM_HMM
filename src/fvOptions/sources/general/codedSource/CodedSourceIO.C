/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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
    applied_.resize(fieldNames_.size(), false);

    dict.readCompat<word>("name", {{"redirectType", 1706}}, name_);


    // Code context chunks

    auto& ctx = codedBase::codeContext();

    ctx.readEntry("codeCorrect", codeCorrect_);
    ctx.readEntry("codeAddSup", codeAddSup_);

    // ctx.readEntry("codeConstrain", codeConstrain_);
    ctx.readEntry  // Compatibility
    (
        coeffs_.lookupEntryCompat
        (
            "codeConstrain",
            {{ "codeSetValue", 1812 }},
            keyType::LITERAL
        ).keyword(),
        codeConstrain_
    );

    return true;
}


// ************************************************************************* //
