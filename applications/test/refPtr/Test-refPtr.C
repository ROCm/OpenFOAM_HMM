/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Application
    Test-refPtr

Description
    Tests some basic functionality of refPtr

\*---------------------------------------------------------------------------*/

#include "primitiveFields.H"
#include "Switch.H"

using namespace Foam;

struct myScalarField : public scalarField
{
    using scalarField::scalarField;
};


template<class T>
void printInfo(const refPtr<T>& item, const bool verbose = false)
{
    Info<< "refPtr good:" << Switch::name(item.good())
        << " pointer:" << Switch::name(item.is_pointer())
        << " addr: " << Foam::name(item.get())
        << " movable:" << Switch(item.movable());

    Info<< " move-constructible:"
        << std::is_move_constructible<refPtr<T>>::value
        << " move-assignable:"
        << std::is_move_assignable<refPtr<T>>::value
        << " nothrow:"
        << std::is_nothrow_move_assignable<refPtr<T>>::value
        << " trivially:"
        << std::is_trivially_move_assignable<refPtr<T>>::value
        << nl;

    if (verbose && item)
    {
        Info<< "content: " << item() << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    {
        Info<< nl << "Construct from reference" << nl;

        scalarField f2(10, Foam::sqrt(2.0));
        printInfo(refPtr<scalarField>(f2), true);
    }

    {
        Info<< nl << "Construct from New (is_pointer)" << nl;
        auto tfld1 = refPtr<scalarField>::New(10, scalar(1));
        printInfo(tfld1, true);

        Info<< nl << "Dereferenced: " << *tfld1 << nl;

        Info<< nl << "Construct from autoPtr" << nl;
        refPtr<scalarField> tfld2(autoPtr<scalarField>::New(10, scalar(2)));
        printInfo(tfld2, true);


        Info<< nl << "Construct from unique_ptr" << nl;
        std::unique_ptr<scalarField> ptr(new scalarField(10, scalar(3)));
        refPtr<scalarField> tfld3(std::move(ptr));
        printInfo(tfld3, true);


        Info<< nl << "Reset from autoPtr" << nl;
        tfld2.reset(autoPtr<scalarField>::New(3, scalar(13)));
        printInfo(tfld2, true);


        Info<< nl << "Reset from unique_ptr" << nl;
        ptr.reset(new scalarField(5, scalar(15)));
        tfld3.reset(std::move(ptr));
        printInfo(tfld3, true);


        ptr.reset(new scalarField(2, scalar(1)));
        Info<< nl << "const-ref from pointer: " << name(ptr.get()) << nl;
        tfld3.cref(ptr.get());
        printInfo(tfld3, true);
    }

    Info<< "\nEnd" << endl;

    return 0;
}


// ************************************************************************* //
