/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Test behaviour of UPtrList, PtrList
\*---------------------------------------------------------------------------*/

#include "OSspecific.H"

#include "scalar.H"
#include "IOstreams.H"
#include "PtrList.H"
#include "DLPtrList.H"
#include "SLPtrList.H"
#include "plane.H"
#include "DynamicList.H"

using namespace Foam;

class Scalar
{
    scalar data_;

public:

    Scalar()
    :
        data_(0)
    {}

    Scalar(scalar val)
    :
        data_(val)
    {}

    ~Scalar()
    {
        Info<<"delete Scalar: " << data_ << endl;
    }

    autoPtr<Scalar> clone() const
    {
        return autoPtr<Scalar>::New(data_);
    }

    friend Ostream& operator<<(Ostream& os, const Scalar& val)
    {
        os  << val.data_;
        return os;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    PtrList<Scalar> list1(10);
    PtrList<Scalar> list2(15);
    PtrList<Scalar> listApp;

    {
        DLPtrList<Scalar> llist1;
        llist1.insert(new Scalar(100));
        llist1.insert(new Scalar(200));
        llist1.insert(new Scalar(300));

        DLPtrList<Scalar>::const_iterator citer = llist1.begin();

        Info<< *citer << endl;
        Info<< typeid(*citer).name() << endl;

        ++citer;
        ++citer;

        --citer;

        Info<< typeid(llist1.begin()).name() << endl;

        forAllIters(llist1, it)
        {
            Info<< typeid(*it).name() << nl
                << "reversed: " << *it << endl;
        }
        for (const auto& it : llist1)
        {
            Info<< typeid(it).name() << nl
                << "for-: " << it << endl;
        }
    }

    forAll(list1, i)
    {
        list1.set(i, new Scalar(1.3*i));
    }

    forAll(list2, i)
    {
        list2.set(i, new Scalar(10 + 1.3*i));
    }

    for (label i = 0; i < 5; ++i)
    {
        listApp.append(new Scalar(1.3*i));
    }

    Info<< nl
        <<"list1: " << list1 << nl
        <<"list2: " << list2 << nl
        <<"list-appended: " << listApp << endl;

    Info<<"indirectly delete some items via set(.., 0) :" << endl;
    for (label i = 0; i < 3; i++)
    {
        list1.set(i, 0);
    }

    Info<<"transfer list2 -> list1:" << endl;
    list1.transfer(list2);

    Info<<"list1: " << list1 << nl
        <<"list2: " << list2 << endl;

    Info<<"indirectly delete some items via setSize :" << endl;
    list1.setSize(4);

    Info<<"list1: " << list1 << endl;

    PtrList<Scalar> list3(std::move(list1));
    Info<<"Move constructed" << endl;

    Info<<"list1: " << list1 << nl
        <<"list2: " << list2 << nl
        <<"list3: " << list3 << endl;


    Info<<"Move construct:" << endl;

    PtrList<Scalar> list4(std::move(list3));

    Info<<"list3: " << list3 << nl
        <<"list4: " << list4 << endl;

    Info<<"Move assign:" << endl;
    list3 = std::move(list4);

    Info<<"list3: " << list3 << nl
        <<"list4: " << list4 << endl;


    Info<<"UPtrList from PtrList" << nl;

    UPtrList<Scalar> ulist1(list3);

    Info<<"ulist1: " << ulist1 << nl;

    Info<<"Move construct:" << endl;

    UPtrList<Scalar> ulist2(std::move(ulist1));

    Info<<"ulist1: " << ulist1 << nl
        <<"ulist2: " << ulist2 << nl;

    Info<<"Copy assign:" << endl;
    ulist1 = ulist2;

    Info<<"ulist1: " << ulist1 << nl
        <<"ulist2: " << ulist2 << nl;

    Info<<"Move assign:" << endl;
    ulist1 = std::move(ulist2);

    Info<<"ulist1: " << ulist1 << nl
        <<"ulist2: " << ulist2 << nl;

    // Test iterator random access
    {
        auto iter1 = ulist1.begin();
        auto iter2 = iter1 + 3;

        Info<<"begin:" << *iter1 << " (+3):" << *iter2 << nl;
        Info<< "diff= " << (iter1 - iter2) << nl;
        Info<< "iter[2]=" << iter1[2] << nl;
        Info<< "iter1 < iter2 : " << (iter1 < iter2) << nl;
        Info<< "iter1 >= iter2 : " << (iter1 >= iter2) << nl;
    }

    PtrList<plane> planes;
    planes.append(new plane(vector::one, vector::one));
    planes.append(new plane(vector(1,2,3), vector::one));

    Info<< nl << "appended values" << nl;
    for (const plane& p : planes)
    {
        Info<< "    plane " << p << endl;
    }

    Info<< nl << "Done." << endl;
    return 0;
}


// ************************************************************************* //
