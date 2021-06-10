/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018
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

#include "objectRegistry.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class GeoField>
Foam::tmp<GeoField>
Foam::expressions::exprDriver::cfindFieldObject
(
    const objectRegistry& obr,
    const word& fldName
)
{
    tmp<GeoField> tfld;

    tfld.cref(obr.cfindObject<GeoField>(fldName));

    return tfld;
}


template<class TableType>
bool Foam::expressions::exprDriver::readTable
(
    const word& name,
    const dictionary& dict,
    HashTable<TableType>& tbl,
    bool clear
)
{
    if (clear)
    {
        tbl.clear();
    }

    if (!dict.found(name))
    {
        return false;
    }

    ITstream& is = dict.lookup(name);
    List<dictionary> input(is);

    for (const dictionary& d : input)
    {
        tbl.insert(dict.get<word>("name"), TableType(d));
    }

    return true;
}


template<class TableType>
void Foam::expressions::exprDriver::writeTable
(
    Ostream& os,
    const word& name,
    const HashTable<TableType>& tbl
)
{
    if (tbl.size())
    {
        os.writeKeyword(name);
        os  << token::BEGIN_LIST << nl;

        forAllConstIters(tbl, iter)
        {
            os.beginBlock();
            os.writeEntry("name", iter.key());
            (*iter).write(os);
            os.endBlock();
        }
        os  << token::END_LIST
            << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Type Foam::expressions::exprDriver::exprDriver::weightedAverage
(
    const scalarField& wfield,
    const Field<Type>& fld
)
{
    if (isNull(wfield))
    {
        const label n = returnReduce(fld.size(), sumOp<label>());

        // stabilize
        if (!n)
        {
            return Zero;
        }

        return gSum(fld) / scalar(n);
    }

    // #ifdef FULLDEBUG
    // checkSize(wfield, fld);
    // #endif

    const scalar s = gSum(wfield);

    // stabilize
    if (mag(s) < ROOTVSMALL)
    {
        return Zero;
    }

    return gSum(wfield*fld) / s;
}


template<class Type>
Type Foam::expressions::exprDriver::exprDriver::weightedSum
(
    const scalarField& wfield,
    const Field<Type>& fld
)
{
    if (isNull(wfield))
    {
        return gSum(fld);
    }

    // #ifdef FULLDEBUG
    // checkSize(wfield, fld);
    // #endif

    return gSum(wfield*fld);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::exprDriver::getResult(bool wantPointData)
{
    if (!result_.isPointData(wantPointData))
    {
        FatalErrorInFunction
            << "Expected a" << (wantPointData ? " point" : "")
            << " field,  but found a" << (!wantPointData ? " point" : "")
            << " field" << nl
            << exit(FatalError);
    }

    return result_.getResult<Type>();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::expressions::exprDriver::isLocalVariable
(
    const word& name,
    bool wantPointData,
    label expectedSize
) const
{
    DebugInfo
        << "Looking for local" << (wantPointData ? " point" : "")
        << " field name:" << name << " type:"
        << pTraits<Type>::typeName << " size:" << expectedSize;


    bool good = hasVariable(name);

    if (good)
    {
        const exprResult& var = variable(name);

        DebugInfo
            << " - found (" << var.valueType() << ' '
            << var.isPointData() << ')';


        good = (var.isType<Type>() && var.isPointData(wantPointData));

        // Do size checking if requested
        if (good && expectedSize >= 0)
        {
            good = (var.size() == expectedSize);
            reduce(good, andOp<bool>());

            if (debug && !good)
            {
                Info<< " size is";
            }
        }
    }

    DebugInfo << (good ? " good" : " bad") << endl;

    return good;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::exprDriver::newField
(
    const Type& val
) const
{
    return tmp<Field<Type>>::New(size(), val);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::exprDriver::newPointField
(
    const Type& val
) const
{
    return tmp<Field<Type>>::New(pointSize(), val);
}


// ************************************************************************* //
