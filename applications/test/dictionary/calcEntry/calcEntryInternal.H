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

Namespace
    Foam::functionEntries::calcEntryInternal

Description
    This dictionary function entry may or may not do anything particularly
    useful - depending upon what is currently being used to test.

SourceFiles
    calcEntryInternal.C

\*---------------------------------------------------------------------------*/

#ifndef calcEntryInternal_H
#define calcEntryInternal_H

#include "functionEntry.H"
#include "memberFunctionSelectionTables.H"
#include "staticMemberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
namespace calcEntryInternal
{

/*---------------------------------------------------------------------------*\
                     Class calcEntryFunctions Declaration
\*---------------------------------------------------------------------------*/

class scalarFunctions
{
public:

    // Member Function Selectors

        declareStaticMemberFunctionSelectionTable
        (
            scalar,
            scalarFunctions,
            dispatch,
            ParamList,
            (
                const UList<scalar>& param
            ),
            (param)
        );

        //- Calculate
        static scalar dispatch(const word&, const UList<scalar>&);


    static scalar pi_0(const UList<scalar>&);
    static scalar degToRad_1(const UList<scalar>&);
    static scalar radToDeg_1(const UList<scalar>&);
    static scalar cos_1(const UList<scalar>&);
    static scalar sin_1(const UList<scalar>&);
    static scalar pow_2(const UList<scalar>&);
    static scalar log_1(const UList<scalar>&);
    static scalar log10_1(const UList<scalar>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace calcEntryInternal
} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
