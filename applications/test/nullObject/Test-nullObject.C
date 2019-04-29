/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014 OpenFOAM Foundation
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
    Test-nullObject

Description
    Tests of nullObject

\*---------------------------------------------------------------------------*/

#include "nullObject.H"
#include "List.H"
#include "HashSet.H"
#include "faceList.H"
#include "pointField.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

class SimpleClass
{
public:

    //- Null constructor
    SimpleClass()
    {}
};


template<class T>
void printInfo(const UList<T>& list)
{
    std::cout
        << nl
        << "List : addr: " << uintptr_t(&list)
        << " (null: " << isNull(list) << ")" << nl
        << "    size: " << list.size() << " empty: " << list.empty() << nl
        << "    data: " << uintptr_t(list.cdata())
        << " begin=" << uintptr_t(list.begin())
        << " end=" << uintptr_t(list.end()) << nl;

    Info<< list << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{
    // Test pointer and reference to a class

    SimpleClass* ptrToClass = new SimpleClass;
    SimpleClass& refToClass(*ptrToClass);

    std::cout
        << "nullObject addr=" << uintptr_t(&(nullObjectPtr)) << nl
        << "  sizeof(nullObject) = " << sizeof(NullObject::nullObject) << nl
        << "  sizeof(void*) = " << sizeof(void*) << nl
        << "  sizeof(labelList) = " << sizeof(labelList) << nl
        << "  sizeof(wordHashSet) = " << sizeof(wordHashSet) << nl << nl;

    std::cout
        << "nullObject" << nl
        << "  pointer:" << uintptr_t(nullObjectPtr->pointer()) << nl
        << "  value:"   << nullObjectPtr->value() << nl << nl;

    if (notNull(ptrToClass))
    {
        Info<< "Pass: ptrToClass is not null" << nl;
    }
    else
    {
        Info<< "FAIL: refToClass is null" << nl;
    }

    if (notNull(refToClass))
    {
        Info<< "Pass: refToClass is not null" << nl;
    }
    else
    {
        Info<< "FAIL: refToClass is null" << nl;
    }


    // Test pointer and reference to the nullObject

    const SimpleClass* ptrToNull(NullObjectPtr<SimpleClass>());
    const SimpleClass& refToNull(*ptrToNull);

    if (isNull(ptrToNull))
    {
        Info<< "Pass: ptrToNull is null" << nl;
    }
    else
    {
        Info<< "FAIL: ptrToNull is not null" << nl;
    }

    if (isNull(refToNull))
    {
        Info<< "Pass: refToNull is null" << nl;
    }
    else
    {
        Info<< "FAIL: refToNull is not null" << nl;
    }

    // Clean-up
    delete ptrToClass;


    // Test List casting
    {
        labelList list1;
        labelList list2({1, 2, 3});

        printInfo(list1);
        printInfo(list2);
        printInfo(labelList::null());

        printInfo(faceList::null());
        printInfo(pointField::null());
    }

    Info<< nl;

    return 0;
}

// ************************************************************************* //
