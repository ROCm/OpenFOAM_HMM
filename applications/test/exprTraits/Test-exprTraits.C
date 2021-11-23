/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Basic tests of expression traits

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "ITstream.H"
#include "exprTraits.H"
#include "uLabel.H"
#include "error.H"
#include "stringList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
void printTraits()
{
    const auto typeCode = exprTypeTraits<T>::value;

    Info<< "type " << pTraits<T>::typeName
        << " code:" << int(typeCode)
        << " name:" << exprTypeTraits<T>::name;

    if (pTraits<T>::typeName != word(exprTypeTraits<T>::name))
    {
        Info<< " (UNSUPPORTED)";
    }

    Info << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{
    Info<< nl << "Traits:" << nl;

    printTraits<word>();
    printTraits<string>();
    printTraits<bool>();
    printTraits<label>();
    printTraits<scalar>();
    printTraits<vector>();
    printTraits<tensor>();
    printTraits<symmTensor>();
    printTraits<sphericalTensor>();

    const auto getName = nameOp<expressions::valueTypeCode>();

    Info<< nl;

    Info<< "Name of typeCode: "
        << Foam::name(expressions::valueTypeCode::type_scalar) << nl;

    Info<< "Name of typeCode: "
        << getName(expressions::valueTypeCode::type_bool) << nl;


    Info<< nl << "Done" << nl;
    return 0;
}


// ************************************************************************* //
