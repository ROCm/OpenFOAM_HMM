/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

Description
    Functionality of IjkField

\*---------------------------------------------------------------------------*/

#include "point.H"
#include "IjkField.H"
#include "IOstreams.H"

using namespace Foam;

template<class T>
Ostream& print(const IjkField<T>& fld)
{
    Info<< static_cast<const Field<T>&>(fld).size()
        << " addr:" << long(fld.cdata()) << ' ' << fld.sizes() << ' '
        << flatOutput(fld);

    return Info;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Create with inconsistent sizes
    IjkField<label> field0({1, 2, 3}, identity(10));

    IjkField<scalar> field1({3, 4, 5});
    IjkField<scalar> field2({2, 3, 3});

    forAll(field1, i)
    {
        field1[i] = -i;
    }
    forAll(field2, i)
    {
        field2[i] = i;
    }

    Info<< "ijk field "; print(field1) << nl;
    Info<< "ijk field "; print(field2) << nl;

    field1.resize(field2.sizes());

    Info<< "resized "; print(field1) << nl;

    field1 *= 2;

    Info<< "Multiply: "; print(field1) << nl;

    field1.resize({1, 2, 3});

    Info<< "Resize - shrink: "; print(field1) << nl;

    field1.resize({2, 3, 2});

    Info<< "Resize - grow: "; print(field1) << nl;

    field1.resize({3, 2, 2});

    Info<< "Resize - repartition: "; print(field1) << nl;

    field1 = field2;

    Info<< "Copied: "; print(field1) << nl;

    field1 = 3.14159;

    Info<< "Assigned: "; print(field1) << nl;

    field1 += 3.14159;

    Info<< "+= operator: "; print(field1) << nl;

    field1 /= 1.2;

    Info<< "/= operator: "; print(field1) << nl;

    IjkField<scalar> field3(std::move(field2));

    Info<< "Move construct: "; print(field2) << nl;
    print(field3) << nl;

    // Field operations are still limited, but we can bypass things too

    {
        Field<scalar>& tmpField = field1;
        tmpField = sqr(tmpField);

        Info<< "squared (workaround): "; print(field1) << nl;
    }


    Info<< nl
        << "Before transfer: addr:" << long(field1.data())
        << " size:" << field1.size() << nl;

    Field<scalar> sfield1(std::move(field1));
    field1.clear();

    Info<< "After transfer to regular field" << nl
        << "    source:" << long(field1.data()) << nl
        << "    target:" << long(sfield1.data()) << nl
        << "Values"
        << "    source:";
    print(field1) << nl;

    Info<< "    target:" << flatOutput(sfield1) << nl;


    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
