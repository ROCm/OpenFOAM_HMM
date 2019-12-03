/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider <bgschaid@hfd-research.com>
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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
Foam::expressions::exprDriver::getResult(bool isPointVal)
{
    if (!result_.isPointValue(isPointVal))
    {
        FatalErrorInFunction
            << "Expected a" << (isPointVal ? " point" : "")
            << " field,  but found a" << (!isPointVal ? " point" : "")
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
    bool isPointVal,
    label expectedSize
) const
{
    DebugInfo
        << "Looking for local" << (isPointVal ? " point" : "")
        << " field name:" << name << " type:"
        << pTraits<Type>::typeName << " size:" << expectedSize;

    bool good = hasVariable(name);

    if (good)
    {
        const exprResult& var = variable(name);

        DebugInfo
            << " - found (" << var.valueType() << ' ' << var.isPointValue() << ')';


        good = (var.isType<Type>() && var.isPointValue(isPointVal));

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
