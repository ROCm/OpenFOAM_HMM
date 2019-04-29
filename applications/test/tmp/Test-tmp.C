/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    Test-tmp

Description
    Tests for possible memory leaks in the tmp (and tmp<Field> algebra).

\*---------------------------------------------------------------------------*/

#include "primitiveFields.H"

using namespace Foam;

struct myScalarField : public scalarField
{
    using scalarField::scalarField;
};


template<class T>
void printInfo(const tmp<T>& tmpItem)
{
    Info<< "tmp valid:" << tmpItem.valid()
        << " isTmp:" << tmpItem.isTmp()
        << " addr: " << uintptr_t(tmpItem.get());

    if (tmpItem.valid())
    {
        Info<< " refCount:" << tmpItem->count();
    }

    Info<< nl;
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

        Info<<"f1 = " << f1 << nl;
    }

    {
        auto tfld1 = tmp<scalarField>::New(20, Zero);

        printInfo(tfld1);

        if (tfld1.valid())
        {
            Info<<"tmp: " << tfld1() << nl;
        }

        // Hold on to the old content for a bit

        tmp<scalarField> tfld2 =
            tmp<scalarField>::NewFrom<myScalarField>(20, Zero);

        printInfo(tfld2);
        if (tfld2.valid())
        {
            Info<<"tmp: " << tfld2() << nl;
        }

        tfld2.clear();

        Info<<"After clear : ";
        printInfo(tfld2);

        tfld2.cref(f1);

        Info<<"Reset const-ref : ";
        printInfo(tfld2);
    }

    Info<< "\nEnd" << endl;
}


// ************************************************************************* //
