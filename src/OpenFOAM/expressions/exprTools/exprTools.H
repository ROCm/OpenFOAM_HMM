/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

Namespace
    Foam::exprTools

Description
    Helper methods for creating expressions.

SourceFiles
    exprTools.C
    exprToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef exprTools_H
#define exprTools_H

#include "exprString.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Namespace exprTools Declaration
\*---------------------------------------------------------------------------*/

namespace exprTools
{

//- Stringified value for zero of given type
template<class Type>
string zeroValue();

//- Stringified version of data with comma separators between components
template<class Type>
string toString(const Type& data, const word& prefix);

//- Stringified version of data with comma separators between components.
//- Uses prefix corresponding to the pTraits of Type.
template<class Type>
string toString(const Type& data);

//- Stringified version of primitiveEntry of Type.
//- Uses comma separators and prefix corresponding to the pTraits of Type
template<class Type>
string toString(ITstream& is);

//- Get an expression string list from a dictionary
//  Supports list and inline specifications, as well as indirect
//  inclusion of other expression string lists.
//
//  Syntax examples,
//  \code
//
//      // Single items. Inline or list format.
//      list1a  "var1=p";
//      list1b  ( "var1=p" );
//
//      // Multiple items. Inline or list format.
//      list2a  "var2=T; var3=U; var4=rho";
//      list2b  ( "var2=T" "var3=U; var4=rho" );
//
//      // Indirection. Inline or list format.
//      list3a  "#list1a; #list2b";
//      list2b  ( "var1=p" "#list2a" );
//
//  \endcode
List<expressions::exprString> getList
(
    const dictionary& dict,
    const word& keyword,
    bool mandatory = true,
    label recursionDepth = 0
);


// Boolean-like values

//- Convert [0-1] values (usually scalars) as false/true
template<class T>
struct toBoolOp
{
    int operator()(const T& val) const
    {
        return (0.5 < Foam::mag(val));
    }
};

//- No conversion needed for bool
template<>
struct toBoolOp<bool>
{
    int operator()(bool val) const
    {
        return val;
    }
};

//- Inverse test of toBoolOp()
template<class T>
struct toBoolNotOp
{
    int operator()(const T& val) const
    {
        return (Foam::mag(val) < 0.5);
    }
};

//- Inverse test of toBoolOp, no conversion needed for bool
template<>
struct toBoolNotOp<bool>
{
    int operator()(bool val) const
    {
        return !val;
    }
};

//- Logical 'and' with scalars as (false, true) values
template<class T>
struct toBoolAndOp
{
    int operator()(const T& a, const T& b) const
    {
        return (toBoolOp<T>()(a) && toBoolOp<T>()(b));
    }
};

//- Logical 'or' with scalars as (false, true) values
template<class T>
struct toBoolOrOp
{
    int operator()(const T& a, const T& b) const
    {
        return (toBoolOp<T>()(a) || toBoolOp<T>()(b));
    }
};

//- Logical 'xor' with scalars as (false, true) values
template<class T>
struct toBoolXorOp
{
    int operator()(const T& a, const T& b) const
    {
        return (toBoolOp<T>()(a) ? toBoolNotOp<T>()(b) : toBoolOp<T>()(b));
    }
};


} // End namespace exprTools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "exprToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
