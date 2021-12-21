/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Application
    Test-tmp

Description
    Tests for possible memory leaks in the tmp (and tmp<Field> algebra).

\*---------------------------------------------------------------------------*/

#include "primitiveFields.H"
#include "Switch.H"

using namespace Foam;

struct myScalarField : public scalarField
{
    using scalarField::scalarField;
};


template<class T>
void printInfo(const tmp<T>& item, const bool verbose = false)
{
    Info<< "tmp good:" << Switch::name(item.good())
        << " pointer:" << Switch::name(item.is_pointer())
        << " addr: " << Foam::name(item.get())
        << " movable:" << Switch(item.movable());
    if (item)
    {
        Info<< " refCount:" << item->count();
    }

    Info<< " move-constructible:"
        << std::is_move_constructible<tmp<T>>::value
        << " move-assignable:"
        << std::is_move_assignable<tmp<T>>::value
        << " nothrow:"
        << std::is_nothrow_move_assignable<tmp<T>>::value
        << " trivially:"
        << std::is_trivially_move_assignable<tmp<T>>::value
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
    scalarField f1(1000000, 1.0), f2(1000000, 2.0), f3(1000000, 3.0);

    {
        for (int iter=0; iter < 50; ++iter)
        {
            f1 = f2 + f3 + f2 + f3;
        }

        Info<< "f1 = " << f1 << nl;
    }

    {
        auto tfld1 = tmp<scalarField>::New(20, Zero);
        printInfo(tfld1, true);

        // Hold on to the old content for a bit

        tmp<scalarField> tfld2 =
            tmp<scalarField>::NewFrom<myScalarField>(20, Zero);

        printInfo(tfld2, true);

        tfld2.clear();

        Info<< "After clear : ";
        printInfo(tfld2);

        tfld2.cref(f1);

        Info<< "Reset to const-ref : ";
        printInfo(tfld2);

        Info<< "Clear const-ref does not affect tmp: ";
        tfld2.clear();
        printInfo(tfld2);

        Info<< "Reset const-ref affects tmp: ";
        tfld2.reset();
        printInfo(tfld2);


        // Reset tfld2 from tfld1
        Info<< "Fld2 : ";
        printInfo(tfld2);

        for (label i = 0; i < 2; ++i)
        {
            tfld2.reset
            (
                tmp<scalarField>::NewFrom<myScalarField>(4, Zero)
            );

            Info<< "Reset to some other tmp content : ";
            printInfo(tfld2);
        }

        std::unique_ptr<scalarField> ptr(new scalarField(2, scalar(1)));

        Info<< nl << "const-ref from pointer: " << name(ptr.get()) << nl;
        tfld2.cref(ptr.get());
        printInfo(tfld2);
    }

    Info<< "\nEnd" << endl;

    return 0;
}


// ************************************************************************* //
