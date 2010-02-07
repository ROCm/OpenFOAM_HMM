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
#include "addToGlobalFunctionSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
namespace calcEntryInternal
{

defineGlobalFunctionSelectionTable(dispatch,ParamList);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define globalConstant0(Name, Constant)\
scalar Name##_0(const UList<scalar>& param)                                   \
{                                                                             \
    return Constant;                                                          \
}                                                                             \
addNamedToGlobalFunctionSelectionTable(dispatch,ParamList,Name##_0,&Name##_0)


#define globalFunction0(Name, Function)\
scalar Name##_0(const UList<scalar>& param)                                   \
{                                                                             \
    return Function();                                                        \
}                                                                             \
addNamedToGlobalFunctionSelectionTable(dispatch,ParamList,Name##_0,&Name##_0)


#define globalFunction1(Name, Function)\
scalar Name##_1(const UList<scalar>& param)                                   \
{                                                                             \
    return Function(param[0]);                                                \
}                                                                             \
addNamedToGlobalFunctionSelectionTable(dispatch,ParamList,Name##_1,&Name##_1)


#define globalFunction2(Name, Function)\
scalar Name##_2(const UList<scalar>& param)                                   \
{                                                                             \
    return Function(param[0], param[1]);                                      \
}                                                                             \
addNamedToGlobalFunctionSelectionTable(dispatch,ParamList,Name##_2,&Name##_2)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

globalConstant0(pi, constant::mathematical::pi);

globalFunction1(degToRad, degToRad);
globalFunction1(radToDeg, radToDeg);
globalFunction1(asin, Foam::asin);
globalFunction1(acos, Foam::acos);
globalFunction1(atan, Foam::atan);
globalFunction1(sin, Foam::sin);
globalFunction1(cos, Foam::cos);
globalFunction1(tan, Foam::tan);
globalFunction1(log, Foam::log);
globalFunction1(log10, Foam::log10);
globalFunction1(mag, Foam::mag);

globalFunction2(atan2, Foam::atan2);
globalFunction2(pow, Foam::pow);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar dispatch(const word& name, const UList<scalar>& param)
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace calcEntryInternal
} // End namespace functionEntries
} // End namespace Foam

// ************************************************************************* //
