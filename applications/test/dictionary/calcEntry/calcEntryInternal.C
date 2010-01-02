/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "calcEntryInternal.H"
#include "addToMemberFunctionSelectionTable.H"
#include "addToStaticMemberFunctionSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
namespace calcEntryInternal
{

defineStaticMemberFunctionSelectionTable(scalarFunctions,dispatch,ParamList);


scalar scalarFunctions::dispatch
(
    const word& name,
    const UList<scalar>& param
)
{
    // create lookup name with parameter count
    const word lookupName = name + '_' + Foam::name(param.size());

    dispatchParamListMemberFunctionTable::iterator mfIter =
        dispatchParamListMemberFunctionTablePtr_->find(lookupName);

    if (mfIter == dispatchParamListMemberFunctionTablePtr_->end())
    {
        FatalErrorIn
        (
            "calcEntryInternal::scalarFunctions::dispatch"
            "(const word&, const UList<scalar>&) : "
        )   << "Unknown function " << name << nl << nl
            << "Valid types are :" << endl
            << dispatchParamListMemberFunctionTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return mfIter()(param);
}


scalar scalarFunctions::pi_0(const UList<scalar>& param)
{
    return constant::mathematical::pi;
}

scalar scalarFunctions::degToRad_1(const UList<scalar>& param)
{
    return degToRad(param[0]);
}

scalar scalarFunctions::radToDeg_1(const UList<scalar>& param)
{
    return radToDeg(param[0]);
}

scalar scalarFunctions::sin_1(const UList<scalar>& param)
{
    return Foam::sin(param[0]);
}

scalar scalarFunctions::cos_1(const UList<scalar>& param)
{
    return Foam::cos(param[0]);
}

scalar scalarFunctions::pow_2(const UList<scalar>& param)
{
    return Foam::pow(param[0], param[1]);
}

scalar scalarFunctions::log_1(const UList<scalar>& param)
{
    return Foam::log(param[0]);
}

scalar scalarFunctions::log10_1(const UList<scalar>& param)
{
    return Foam::log10(param[0]);
}



addNamedToStaticMemberFunctionSelectionTable
(
    scalarFunctions,scalarFunctions,dispatch,ParamList,
    pi_0,
    &scalarFunctions::pi_0
);

addNamedToStaticMemberFunctionSelectionTable
(
    scalarFunctions,scalarFunctions,dispatch,ParamList,
    degToRad_1,
    &scalarFunctions::degToRad_1
);

addNamedToStaticMemberFunctionSelectionTable
(
    scalarFunctions,scalarFunctions,dispatch,ParamList,
    radToDeg_1,
    &scalarFunctions::radToDeg_1
);

addNamedToStaticMemberFunctionSelectionTable
(
    scalarFunctions,scalarFunctions,dispatch,ParamList,
    sin_1,
    &scalarFunctions::sin_1
);

addNamedToStaticMemberFunctionSelectionTable
(
    scalarFunctions,scalarFunctions,dispatch,ParamList,
    cos_1,
    &scalarFunctions::cos_1
);

addNamedToStaticMemberFunctionSelectionTable
(
    scalarFunctions,scalarFunctions,dispatch,ParamList,
    pow_2,
    &scalarFunctions::pow_2
);

addNamedToStaticMemberFunctionSelectionTable
(
    scalarFunctions,scalarFunctions,dispatch,ParamList,
    log_1,
    &scalarFunctions::log_1
);



addNamedToStaticMemberFunctionSelectionTable
(
    scalarFunctions,scalarFunctions,dispatch,ParamList,
    log10_1,
    &scalarFunctions::log10_1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace calcEntryInternal
} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
