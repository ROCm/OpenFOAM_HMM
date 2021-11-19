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
#include "exprScanToken.H"

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


    {
        expressions::scanToken tok;
        expressions::scanToken tok2;

        Info<< nl << "sizeof(scanToken): "
            << sizeof(tok) << nl;

        Info<< "    type:" << int(tok.type_) << nl;
        Info<< "    ptr:" << Foam::name(tok.name_) << nl;

        Info<< "    type:" << int(tok2.type_) << nl;
        Info<< "    ptr:" << Foam::name(tok2.name_) << nl;

        tok.setWord("hello");

        Info<< "    type:" << int(tok.type_) << nl;
        Info<< "    ptr:" << Foam::name(tok.name_) << nl;

        tok2 = tok;
        Info<< "    type:" << int(tok2.type_) << nl;
        Info<< "    ptr:" << Foam::name(tok2.name_) << nl;

        tok2.destroy();

        Info<< "    type:" << int(tok2.type_) << nl;
        Info<< "    ptr:" << Foam::name(tok2.name_) << nl;
    }

    Info<< nl << "Done" << nl;
    return 0;
}


// ************************************************************************* //
